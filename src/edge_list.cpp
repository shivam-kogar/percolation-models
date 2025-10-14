
#include <vector>
#include <iostream>
#include "is_edge.hpp"
#include "is_site.hpp"
#include "edge_list.hpp"
#include "graph_utils.hpp"

using namespace std;

bool is_edge(double p);
bool is_site(double p);

vector< pair<int, int> > edge_list_2d_bond(double p, int n) {
    vector<pair<int, int> > edge_list; 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int id_curr = pair_to_id(i, j, n);
            int id_down = pair_to_id(i + 1, j, n);
            int id_right = pair_to_id(i, j + 1, n);
            if (i + 1 < n) {
                if (is_edge(p)) {
                    edge_list.emplace_back(id_curr, id_down);
                }
            }
            if (j + 1 < n) {
                if (is_edge(p)) {
                    edge_list.emplace_back(id_curr, id_right);   
                }
            }
        }
    }
    return edge_list;
}         

vector<pair<int, int>> edge_list_2d_site(double p, int n) {
    vector<pair<int, int>> edge_list;
    vector<vector<bool>> site_open(n, vector<bool>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            site_open[i][j] = is_site(p);  
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!site_open[i][j]) continue;

            int id_curr = pair_to_id(i, j, n);

            if (i + 1 < n && site_open[i + 1][j]) {
                int id_down = pair_to_id(i + 1, j, n);
                edge_list.emplace_back(id_curr, id_down);
            }

            if (j + 1 < n && site_open[i][j + 1]) {
                int id_right = pair_to_id(i, j + 1, n);
                edge_list.emplace_back(id_curr, id_right);
            }
        }
    }

    return edge_list;
}

