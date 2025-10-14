#include <vector>
#include <unordered_set>
#include <utility>
#include <cmath>
#include "graph_utils.hpp"

using namespace std;

struct pair_hash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

bool within_bounds(int ni, int nj, int n) {
    return (0 <= ni && ni < n && 0 <= nj && nj < n);
}

unordered_set<int> ordered_strip_nodes(const vector<int>& geodesic_ids, int dist, int n) {
    unordered_set<int> strip_nodes;
    for (int id : geodesic_ids) {
        auto coord = id_to_pair(id, n);
        for (int dx = -dist; dx <= dist; ++dx) {
            for (int dy = -dist; dy <= dist; ++dy) {
                if (abs(dx) + abs(dy) <= dist) {
                    int ni = coord.first + dx;
                    int nj = coord.second + dy;
                    if (within_bounds(ni, nj, n)) {
                        strip_nodes.insert(pair_to_id(ni, nj, n));
                    }
                }
            }
        }
    }
    return strip_nodes; 
}

vector<pair<int, int>> strip_edges(
    const unordered_set<int>& strip_nodes,
    const vector<pair<int, int>>& edges,
    int n, bool ordered
) {
    vector<pair<int, int>> strip_edges;
    vector<pair<int, int>> directions = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

    unordered_set<pair<int, int>, pair_hash> edge_set;
    for (const auto& e : edges) {
        edge_set.insert(e);
        edge_set.insert({e.second, e.first}); 
    }

    for (int id : strip_nodes) {
        auto coord = id_to_pair(id, n);
        for (auto [di, dj] : directions) {
            int ni = coord.first + di;
            int nj = coord.second + dj;

            if (within_bounds(ni, nj, n)) {
                int nbr_id = pair_to_id(ni, nj, n);
                if (strip_nodes.count(nbr_id) && id < nbr_id) {
                    if (ordered) {
                        strip_edges.emplace_back(id, nbr_id);
                    }
                    else {
                        if (edge_set.count({id, nbr_id})) {
                            strip_edges.emplace_back(id, nbr_id);
                        }
                    }
                }
            }
        }
    }
    return strip_edges;
}

vector<pair<int, int>> get_strip_edges_from_ids(unordered_set<int> strip_nodes,
    const vector<pair<int, int>>& edges, int n, bool ordered
) {
    return strip_edges(strip_nodes, edges, n, ordered);
}

double cross_sectional_width(const vector<int>& geodesic, int n, 
                             const unordered_set<int>& strip_nodes) {
    double width_sum = 0.0;
    int count = 0;

    for (int i = 1; i < geodesic.size() - 1; ++i) {
        auto [x_prev, y_prev] = id_to_pair(geodesic[i - 1], n);
        auto [x_curr, y_curr] = id_to_pair(geodesic[i], n);
        auto [x_next, y_next] = id_to_pair(geodesic[i + 1], n);

        int dx1 = x_curr - x_prev;
        int dy1 = y_curr - y_prev;
        int dx2 = x_next - x_curr;
        int dy2 = y_next - y_curr;

        if (dx1 != 0) dx1 /= abs(dx1);
        if (dy1 != 0) dy1 /= abs(dy1);
        if (dx2 != 0) dx2 /= abs(dx2);
        if (dy2 != 0) dy2 /= abs(dy2);

        auto sweep = [&](int dx_ortho, int dy_ortho) -> int {
            int w = 1;
            for (int s = 1;; ++s) {
                int nx = x_curr + s * dx_ortho;
                int ny = y_curr + s * dy_ortho;
                if (within_bounds(nx, ny, n) && 
                    strip_nodes.count(pair_to_id(nx, ny, n))) {
                    w++;
                } else {
                    break;
                }
            }
            for (int s = 1;; ++s) {
                int nx = x_curr - s * dx_ortho;
                int ny = y_curr - s * dy_ortho;
                if (within_bounds(nx, ny, n) && 
                    strip_nodes.count(pair_to_id(nx, ny, n))) {
                    w++;
                } else {
                    break;
                }
            }
            return w;
        };

        double width;

        if (dx1 != dx2 || dy1 != dy2) {
            // Turning point → average incoming & outgoing orthogonal widths
            int w1 = sweep(-dy1, dx1); // orthogonal to incoming
            int w2 = sweep(-dy2, dx2); // orthogonal to outgoing
            width = 0.5 * (w1 + w2);
        } else {
            // Regular point → use one direction
            width = sweep(-dy1, dx1);
        }

        width_sum += width;
        count++;
    }

    return width_sum / count;
}

int compute_dual_area(const unordered_set<int>& strip_nodes, int n) {
    unordered_set<pair<int, int>, pair_hash> dual_faces;

    for (int id : strip_nodes) {
        auto [i, j] = id_to_pair(id, n);

        for (int di = -1; di <= 0; ++di) {
            for (int dj = -1; dj <= 0; ++dj) {
                int fi = i + di;
                int fj = j + dj;
                if (0 <= fi && fi < n - 1 && 0 <= fj && fj < n - 1) {
                    dual_faces.insert({fi, fj});
                }
            }
        }
    }

    return dual_faces.size(); 
}
