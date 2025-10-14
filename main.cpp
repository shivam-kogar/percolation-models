#include <vector>
#include <sstream>
#include <fstream>
#include <limits>
#include <omp.h>
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


struct Stats {
    double mean;
    double se;
};

struct AllStats {
    Stats geodesic;
    vector<Stats> disordered;
    vector<Stats> ordered;
    Stats exact;
    vector<Stats> area;
    Stats length;
    vector<Stats> width;
    int num_zero;
};

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

struct ProcessedData {
    vector<RawData> data;
    vector<bool> no_path_flags;
    int num_zero;
};

vector<pair<int, int>> edge_list_2d_site(double p, int n);
vector< pair<int, int> > edge_list_2d_bond(double p, int n);
int compute_dual_area(const unordered_set<int>& strip_nodes, int n);
double cross_sectional_width(const vector<int>& geodesic, int n, 
                             const unordered_set<int>& strip_nodes);
unordered_set<int> ordered_strip_nodes(const vector<int>& geodesic_ids, int dist, int n);
vector<pair<int, int>> get_strip_edges_from_ids(unordered_set<int> strip_nodes,
    const vector<pair<int, int>>& edges, int n, bool ordered
);
SpMat laplacian_2d(const vector<pair<int, int> >& edges, int n);
double resistance_distance_2d(int i, int j, const SpMat& L, double s);
ProcessedData simulate_all(int id_start, int id_end, double p);
RawData simulate_conductances(int id_start, int id_end, double p, EdgeGenerator edge_list_fn);
void save_data(string filename, pair<int, int> start, pair<int, int> end);
AllStats make_stats(const ProcessedData& y, bool with_zero);

const int n = 1000; // we are studying diffusion on a square, 2D n x n lattice
const double s = 1e-2; // inverse-time parameter in Laplace Domain
int num_sims = []{
    const char* v = std::getenv("NUM_SIMS");
    return v ? std::max(1, atoi(v)) : 10;  // default to 10
}(); // number of lattice configurations for each value of p
const int max_dist = 5;
const int num_data = 10;
const double p_c = 0.5;
const double ln_min = log(0.5001 - p_c);
const double ln_max = log(0.6 - p_c);

int main(int argc, char** argv) {
    Eigen::setNbThreads(1);
    omp_set_num_threads(omp_get_max_threads());

    int p_index = 0;
    if (argc > 1) p_index = atoi(argv[1]);

    const char* jobid  = getenv("SLURM_JOB_ID");
    const char* taskid = getenv("SLURM_ARRAY_TASK_ID");

    std::ostringstream fname;
    fname << "output_job" << (jobid ? jobid : "local")
          << "_task" << (taskid ? taskid : "0")
          << ".csv";

    save_data(fname.str(), {249,249}, {749,749}, p_index);
}




void save_data(string filename, pair<int, int> start, pair<int, int> end, int p_index) {
    int id_start = pair_to_id(start.first, start.second, n);
    int id_end   = pair_to_id(end.first, end.second, n);

    // compute the p for this index
    double ln_offset = ln_min + p_index * (ln_max - ln_min) / (num_data - 1);
    double p = p_c + exp(ln_offset);
    cout << "Percolation probability: " << p << "\n";

    // run simulation
    const auto y             = simulate_all(id_start, id_end, p);
    const auto y_with_zero   = make_stats(y, true);
    const auto y_without_zero= make_stats(y, false);
    int num_zero             = y_with_zero.num_zero;

    // open CSV for this p
    ofstream out(filename);

    // header (same as your current header code)
    out << "p,geodesic_mean,geodesic_se,";
    for (int j = 0; j < max_dist; ++j) {
        out << "disordered_" << j+1 << "_mean,"
            << "disordered_" << j+1 << "_se,"
            << "ordered_"    << j+1 << "_mean,"
            << "ordered_"    << j+1 << "_se,";
    }
    out << "exact_mean,exact_se,";
    for (int j = 0; j < max_dist; ++j)
        out << "area_" << j+1 << "_mean,area_" << j+1 << "_se,";
    out << "length_mean,length_se,";
    for (int j = 0; j < max_dist; ++j)
        out << "width_" << j+1 << "_mean,width_" << j+1 << "_se,";
    out << "geodesic_mean_zero,geodesic_se_zero,";
    for (int j = 0; j < max_dist; ++j) {
        out << "disordered_" << j+1 << "_mean_zero,"
            << "disordered_" << j+1 << "_se_zero,"
            << "ordered_"    << j+1 << "_mean_zero,"
            << "ordered_"    << j+1 << "_se_zero,";
    }
    out << "exact_mean_zero,exact_se_zero,";
    for (int j = 0; j < max_dist; ++j)
        out << "area_" << j+1 << "_mean_zero,area_" << j+1 << "_se_zero,";
    out << "length_mean_zero,length_se_zero,";
    for (int j = 0; j < max_dist; ++j)
        out << "width_" << j+1 << "_mean_zero,width_" << j+1 << "_se_zero,";
    out << "num_zero\n";

    // data row
    out << p << ","
        << y_without_zero.geodesic.mean << "," << y_without_zero.geodesic.se << ",";
    for (int j = 0; j < max_dist; ++j)
        out << y_without_zero.disordered[j].mean << "," << y_without_zero.disordered[j].se << ","
            << y_without_zero.ordered[j].mean    << "," << y_without_zero.ordered[j].se    << ",";
    out << y_without_zero.exact.mean << "," << y_without_zero.exact.se << ",";
    for (int j = 0; j < max_dist; ++j)
        out << y_without_zero.area[j].mean << "," << y_without_zero.area[j].se << ",";
    out << y_without_zero.length.mean << "," << y_without_zero.length.se << ",";
    for (int j = 0; j < max_dist; ++j)
        out << y_without_zero.width[j].mean << "," << y_without_zero.width[j].se << ",";
    out << y_with_zero.geodesic.mean << "," << y_with_zero.geodesic.se << ",";
    for (int j = 0; j < max_dist; ++j)
        out << y_with_zero.disordered[j].mean << "," << y_with_zero.disordered[j].se << ","
            << y_with_zero.ordered[j].mean    << "," << y_with_zero.ordered[j].se    << ",";
    out << y_with_zero.exact.mean << "," << y_with_zero.exact.se << ",";
    for (int j = 0; j < max_dist; ++j)
        out << y_with_zero.area[j].mean << "," << y_with_zero.area[j].se << ",";
    out << y_with_zero.length.mean << "," << y_with_zero.length.se << ",";
    for (int j = 0; j < max_dist; ++j)
        out << y_with_zero.width[j].mean << "," << y_with_zero.width[j].se << ",";
    out << num_zero << "\n";

    out.close();
    cout << "Saved data to " << filename << "\n";
}


ProcessedData simulate_all(int id_start, int id_end, double p) {
    int num_zero = 0;
    vector<RawData> all_data(num_sims);
    vector<bool> no_path_flags(num_sims, false);  

    #pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < num_sims; i++) {
        uint64_t seed = 0xcbf29ce484222325ULL
                  ^ (uint64_t)i * 0x100000001b3ULL
                  ^ (uint64_t)std::llround((p - 0.5)*1e9);
        auto c = simulate_conductances(id_start, id_end, p, edge_list_2d_bond);
        all_data[i] = c;

        if (c.no_path) {
            no_path_flags[i] = true;
            #pragma omp atomic
                num_zero++;
        }
    }

    return {all_data, no_path_flags, num_zero};
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
/*
pair<double, double> simulate_conductances(int id_start, int id_end, double p) {
    auto edges = edge_list_2d(p, n);
    auto adj = edge_list_to_adjacency_list(edges, n);

    UnionFind uf(n);

    for (const auto& [u, v] : edges) {
        uf.unite(u, v);
    }

    if (!uf.connected(id_start, id_end)) {
        return {0.0, 0.0}; 
    }

    auto path = bfs(id_start, id_end, adj, n);


    sp_mat L = laplacian_2d(edges, n);
    double R = resistance_distance_2d(id_start, id_end, L, s);
    // cout << "\nResistance Distance: " << R << "\n";
    double C = 1 / R;

    double R_g = path.size() - 1;
    double C_g = 1 / R_g;

    return {C, C_g};
}
     */

/** pair<double, double> simulate_conductances(int id_start, int id_end, vector< pair<int, int> > edges) {
    auto adj = edge_list_to_adjacency_list(edges, n);
    auto geodesic = bfs_adaptive(id_start, id_end, adj, n);
    if (geodesic.path.empty()) {
        cerr << "No pathway!";
    }
    auto edges_g = path_to_edge_list(geodesic.path);
    if (id_start != geodesic.actual_start or id_end != geodesic.actual_end) {
        cout << "\n" << "Changed Endpoints from " << id_to_string(id_start, n) << " and " << id_to_string(id_end, n) << 
        " to " << id_to_string(geodesic.actual_start, n) << " and " << id_to_string(geodesic.actual_end, n) << "\n";
    }

    sp_mat L = laplacian_2d(edges, n);
    double R = resistance_distance_2d(geodesic.actual_start, geodesic.actual_end, L, s);
    cout << "\n" << "Resistance Distance: " << R << "\n";
    double C = 1 / R;

    sp_mat L_g = laplacian_2d(edges_g, n);
    double R_g = resistance_distance_2d(geodesic.actual_start, geodesic.actual_end, L_g, s);
    cout << "\n" << "Resistance Distance (Geodesic): " << R_g << "\n";
    double C_g = 1 / R_g;

    return {C, C_g};
} **/

