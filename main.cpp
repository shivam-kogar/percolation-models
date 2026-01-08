#include <vector>
#include <sstream>
#include <fstream>
#include <limits>
#include <iostream>
#include <numeric>
#include <Eigen/Sparse>  
#include "edge_list.hpp"
#include "laplacian.hpp"
#include "graph_utils.hpp"
#include "math_utils.hpp"
#include "strip_utils.hpp"
#include "bfs.hpp"
#include "resistance_distance.hpp"
#include <unordered_set>
#include <functional> 

using namespace std;
using SpMat = Eigen::SparseMatrix<double>;
using EdgeGenerator = function<vector<pair<int, int>>(double, int)>;

struct RawData {
    double geodesic;
    vector<double> disordered;
    vector<double> ordered;
    double exact;
    vector<double> area;
    double length;
    vector<double> width;
    bool no_path;
};

// vector<pair<int, int>> edge_list_2d_site(double p, int n);
vector<pair<int, int>> edge_list_2d_bond(double p, int n); // given bond probability p, generate an edge list that represents a lattice

int compute_dual_area(const unordered_set<int>& strip_nodes, int n); // computes the area of the strip by counting the number of nodes in the strip's dual graph
double cross_sectional_width(const vector<int>& geodesic, int n,
                             const unordered_set<int>& strip_nodes); // computes the width of a strip about the geodesic
unordered_set<int> ordered_strip_nodes(const vector<int>& geodesic_ids, int dist, int n); // generates a list of nodes in the ordered strip about the geodesic, given the geodesic

vector<pair<int, int>> get_strip_edges_from_ids( // creates an edge list pertaining to the strip 
    unordered_set<int> strip_nodes,
    const vector<pair<int, int>>& edges, int n, bool ordered
);

SpMat laplacian_2d(const vector<pair<int, int>>& edges, int n); // generates the Laplacian matrix of the lattice given an edge list representing the lattice
double resistance_distance_2d(int i, int j, const SpMat& L, double s); // computes the resistance between two endpoints across the lattice

RawData simulate_conductances(int id_start, int id_end, double p, EdgeGenerator edge_list_fn); 
vector<RawData> simulate_all(int id_start, int id_end, double p);

void write_raw_csv(const string& filename, double p, const vector<RawData>& sims);

const int n = 1000; // we are studying diffusion on a square, 2D n x n lattice
const double s = 1e-2; // inverse-time parameter in Laplace Domain
int num_sims = 10 // number of lattice configurations for each value of p, adjust as necessary
const int max_dist = 5;
const int num_data = 10; // number of configurations for each value of p
const double p_c = 0.5; // percolation threshold, we are interested in conductance statistics as p -> p_c
const double ln_min = log(0.5001 - p_c);
const double ln_max = log(0.6 - p_c);

int main(int argc, char** argv) {
    int p_index = (argc > 1) ? atoi(argv[1]) : 0;

    // fixed endpoints
    pair<int,int> start{249,249}, end{749,749};
    int id_start = pair_to_id(start.first, start.second, n);
    int id_end   = pair_to_id(end.first, end.second, n);

    // log-spaced p above p_c
    double ln_offset = ln_min + p_index * (ln_max - ln_min) / (num_data - 1);
    double p = p_c + exp(ln_offset);
    cout << "Percolation probability: " << p << "\n";

    // run sims (raw)
    vector<RawData> sims = simulate_all(id_start, id_end, p);

    // write raw CSV (one row per simulation)
    string filename = "raw_pindex_" + to_string(p_index) + ".csv";
    write_raw_csv(filename, p, sims);

    cout << "Saved data to " << filename << "\n";
    return 0;
}

void write_raw_csv(const string& filename, double p, const vector<RawData>& sims) {
    ofstream out(filename);

    // header
    out << "sim,p,no_path,geodesic,exact,length";
    for (int j = 1; j <= max_dist; ++j) out << ",disordered_" << j;
    for (int j = 1; j <= max_dist; ++j) out << ",ordered_" << j;
    for (int j = 1; j <= max_dist; ++j) out << ",area_" << j;
    for (int j = 1; j <= max_dist; ++j) out << ",width_" << j;
    out << "\n";

    // rows
    for (int i = 0; i < (int)sims.size(); ++i) {
        const auto& r = sims[i];
        out << i << "," << p << "," << (r.no_path ? 1 : 0) << ","
            << r.geodesic << "," << r.exact << "," << r.length;

        for (double x : r.disordered) out << "," << x;
        for (double x : r.ordered)    out << "," << x;
        for (double x : r.area)       out << "," << x;
        for (double x : r.width)      out << "," << x;

        out << "\n";
    }
}


vector<RawData> simulate_all(int id_start, int id_end, double p) {
    vector<RawData> all_data;
    all_data.reserve(num_sims);

    for (int i = 0; i < num_sims; ++i) {
        all_data.push_back(simulate_conductances(id_start, id_end, p, edge_list_2d_bond));
    }
    return all_data;
}

AllStats make_stats(const ProcessedData& y, bool with_zero) {
    const auto& sim_data = y.data;
    const auto& flags = y.no_path_flags;

    vector<vector<double>> disordered_vals(max_dist);
    vector<vector<double>> ordered_vals(max_dist);
    vector<vector<double>> area_vals(max_dist);
    vector<vector<double>> width_vals(max_dist);
    vector<double> geodesic_vals;
    vector<double> exact_vals;
    vector<double> length_vals;

    for (size_t i = 0; i < sim_data.size(); ++i) {
        if (!with_zero && flags[i]) continue;

        const RawData& sim = sim_data[i];
        geodesic_vals.push_back(sim.geodesic);
        exact_vals.push_back(sim.exact);
        length_vals.push_back(sim.length);

        for (int j = 0; j < max_dist; ++j) {
            disordered_vals[j].push_back(sim.disordered[j]);
            ordered_vals[j].push_back(sim.ordered[j]);
            area_vals[j].push_back(sim.area[j]);
            width_vals[j].push_back(sim.width[j]);
        }
    }

    int denom = geodesic_vals.size();
    auto compute = [denom](const vector<double>& v) -> Stats {
        if (denom == 0) return {0.0, 0.0};
        double m = mean(v);
        return {m, stdev(v, m) / sqrt(denom)};
    };

    vector<Stats> disordered(max_dist), ordered(max_dist), area(max_dist), width(max_dist);
    for (int i = 0; i < max_dist; ++i) {
        disordered[i] = compute(disordered_vals[i]);
        ordered[i] = compute(ordered_vals[i]);
        area[i] = compute(area_vals[i]);
        width[i] = compute(width_vals[i]);
    }

    Stats geodesic_stat = compute(geodesic_vals);
    Stats exact_stat = compute(exact_vals);
    Stats length_stat = compute(length_vals);

    return {
        geodesic_stat, disordered, ordered,
        exact_stat, area, length_stat, width, y.num_zero
    };
}


RawData simulate_conductances(int id_start, int id_end, double p, EdgeGenerator edge_list_fn) {
    auto edges = edge_list_fn(p, n);
    auto adj = edge_list_to_adjacency_list(edges, n);
    bool no_path = false;

    UnionFind uf(n);

    for (const auto& [u, v] : edges) {
        uf.unite(u, v);
    }

    vector<double> zeroes(max_dist, 0.0);

    if (!uf.connected(id_start, id_end)) {
        no_path = true;
        return {0.0, zeroes, zeroes, 0.0, zeroes, 0.0, zeroes, no_path};
    }

    auto path = bfs(id_start, id_end, adj, n);

    SpMat L_exact = laplacian_2d(edges, n);
    double R_exact = resistance_distance_2d(id_start, id_end, L_exact, s);
    
    double R_geodesic = path.size() - 1;
    double length = R_geodesic;

    vector<unordered_set<int>> ordered_strips(max_dist);
    vector<double> areas(max_dist);
    vector<double> widths(max_dist);
    vector<vector<pair<int, int>>> ordered_strip_edges(max_dist);
    vector<vector<pair<int, int>>> disordered_strip_edges(max_dist);
    vector<SpMat> ordered_laplacians(max_dist);
    vector<SpMat> disordered_laplacians(max_dist);
    vector<double> ordered_conductances(max_dist);
    vector<double> disordered_conductances(max_dist);

    for (int i = 0; i < max_dist; i++) {
        ordered_strips[i] = ordered_strip_nodes(path, i + 1, n);
        areas[i] = compute_dual_area(ordered_strips[i], n);
        widths[i] = cross_sectional_width(path, n, ordered_strips[i]);

        ordered_strip_edges[i] = get_strip_edges_from_ids(ordered_strips[i], edges, n, true);
        disordered_strip_edges[i] = get_strip_edges_from_ids(ordered_strips[i], edges, n, false);

        ordered_laplacians[i] = laplacian_2d(ordered_strip_edges[i], n);
        disordered_laplacians[i] = laplacian_2d(disordered_strip_edges[i], n);

        ordered_conductances[i] = 1.0 / resistance_distance_2d(id_start, id_end, ordered_laplacians[i], s);
        disordered_conductances[i] = 1.0 / resistance_distance_2d(id_start, id_end, disordered_laplacians[i], s);
    }

    return {1.0 / R_geodesic, 
        disordered_conductances, 
        ordered_conductances, 
        1.0 / R_exact, 
        areas, length, widths, no_path};
}
