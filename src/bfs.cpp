#include <queue>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "graph_utils.hpp"

using namespace std;

vector<int> reconstruct_path(int meeting_node, int id_start, int id_end, 
    const vector<int>& parents_start, const vector<int>& parents_end);

/**
struct AdaptedPath {
    vector<int> path;
    int actual_start;
    int actual_end;
};

vector<int> bfs(int id_start, int id_end, const vector<vector<int>>& adj, int n);
vector<int> neighbors(int id, int n);

 AdaptedPath bfs_adaptive(int id_start, int id_end, const vector<vector<int>>& adj, int n) {
    auto path = bfs(id_start, id_end, adj, n);
    if (!path.empty()) return {path, id_start, id_end};

    auto nbrs_start = neighbors(id_start, n);
    for (int nbr : nbrs_start) {
        auto path1 = bfs(nbr, id_end, adj, n);
        if (!path1.empty()) return {path1, nbr, id_end};
    }

    auto nbrs_end = neighbors(id_end, n);
    for (int nbr : nbrs_end) {
        auto path2 = bfs(id_start, nbr, adj, n);
        if (!path2.empty()) return {path2, id_start, nbr};
    }

    for (int nbr_start : nbrs_start) {
        for (int nbr_end : nbrs_end) {
            auto path3 = bfs(nbr_start, nbr_end, adj, n);
            if (!path3.empty()) return {path3, nbr_start, nbr_end};
        }
    }

    cerr << "No Path Exists, Even Between Neighbors!\n";
    return {};
} 

vector<int> neighbors(int id, int n) {
    vector<int> neighbors;
    auto nbr = id_to_pair(id, n);
    if (nbr.first - 1 > -1) {
        neighbors.push_back(pair_to_id(nbr.first - 1, nbr.second, n));
    }
    if (nbr.first + 1 < n) {
        neighbors.push_back(pair_to_id(nbr.first + 1, nbr.second, n));
    }
    if (nbr.second - 1 > -1) {
        neighbors.push_back(pair_to_id(nbr.first, nbr.second - 1, n));
    }
    if (nbr.second + 1 < n) {
        neighbors.push_back(pair_to_id(nbr.first, nbr.second + 1, n));
    }
    return neighbors;
} **/

vector<int> bfs(int id_start, int id_end, const vector<vector<int>>& adj, int n) {
    int N = n * n; // number of nodes after 2D lattice has been flattened
    queue<int> to_visit_start;
    vector<int> parents_start(N, -1);
    vector<bool> visited_start(N, false); 

    queue<int> to_visit_end;
    vector<int> parents_end(N, -1);
    vector<bool> visited_end(N, false); 

    int meeting_node = -1;

    bool found = false;

    to_visit_start.push(id_start);
    to_visit_end.push(id_end);
    visited_start[id_start] = true;
    visited_end[id_end] = true;
    while (!to_visit_start.empty() and !to_visit_end.empty()) {
        int x = to_visit_start.front();
        to_visit_start.pop();
        for (int i : adj[x]) {
            if (!visited_start[i]) {
                to_visit_start.push(i);
                visited_start[i] = true;
                parents_start[i] = x;
                if (visited_end[i]) {
                    meeting_node = i;
                    found = true;
                    break;
                }
            }
        }
        if (found) {
            break;
        }
        int y = to_visit_end.front();
        to_visit_end.pop();
        for (int i : adj[y]) {
            if (!visited_end[i]) {
                to_visit_end.push(i);
                visited_end[i] = true;
                parents_end[i] = y;
                if (visited_start[i]) {
                    meeting_node = i;
                    found = true;
                    break;
                }
            }
        }
        if (found) {
            break;
        }
    }
    return reconstruct_path(meeting_node, id_start, id_end, parents_start, parents_end);
}

vector<int> reconstruct_path(int meeting_node, int id_start, int id_end, 
    const vector<int>& parents_start, const vector<int>& parents_end) {
    vector<int> path;

    int node = meeting_node;
    while (node != id_start) {
        path.push_back(node);
        node = parents_start[node];
    }
    path.push_back(id_start);
    reverse(path.begin(), path.end());

    node = meeting_node;
    while (node != id_end) {
        node = parents_end[node];
        path.push_back(node);
    }

    return path;
}
