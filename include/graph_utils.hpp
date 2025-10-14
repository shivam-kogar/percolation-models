#pragma once
#include <utility>
#include <string>
#include <vector>
#include <numeric>
#include <sstream>
using namespace std;

class UnionFind {
    public:
        vector<int> parent;
        vector<int> rank;

    UnionFind(int n) {
        int N = n * n;
        parent.resize(N);
        rank.resize(N);
        for (int i = 0; i < N; i++) {
            parent[i] = i;
        }
    }

    int find(int u) {
        if (parent[u] != u) {
            parent[u] = find(parent[u]);
        }
        return parent[u];
    }

    void unite(int u, int v) {
        int root_u = find(u);
        int root_v = find(v);
        if (rank[root_u] < rank[root_v]) {
            parent[root_u] = root_v;
        } 
        else if (rank[root_u] > rank[root_v]) {
            parent[root_v] = root_u;
        } 
        else {
            parent[root_v] = root_u;
            rank[root_u]++;
        }
    }

    bool connected(int u, int v) {
        return find(u) == find(v);
    }
};

inline int pair_to_id(int i, int j, int n) {
    return i * n + j;
}

inline pair<int, int> id_to_pair(int id, int n) {
    return {id / n, id % n};
}

inline string pair_to_string(const pair<int, int>& p) {
    return "(" + to_string(p.first) + ", " + to_string(p.second) + ")";
}

inline string id_to_string(int id, int n) {
    return pair_to_string(id_to_pair(id, n));
}

inline vector<pair<int, int>> path_to_edge_list(const vector<int>& path) {
    vector<pair<int, int>> edges;
    if (path.size() < 2) {
        return edges;
    }
    for (int i = 0; i < path.size() - 1; i++) {
        edges.emplace_back(path[i], path[i + 1]);
    }
    return edges;
}

inline vector<vector<int>> edge_list_to_adjacency_list(const vector<pair<int, int>>& edges, int n) {
    int N = n * n;
    vector<vector<int>> adj(N);
    for (auto edge : edges) {
        adj[edge.first].push_back(edge.second);
        adj[edge.second].push_back(edge.first);
    }
    return adj;
}

inline string vector_to_string(const vector<double>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i != vec.size() - 1) oss << ", ";
    }
    return oss.str();
}

inline string vector_pair_to_string(const vector<pair<int, int>>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << "\n (" << vec[i].first << ", " << vec[i].second << ")\n";
        if (i != vec.size() - 1) oss << ", ";
    }
    return oss.str();
}