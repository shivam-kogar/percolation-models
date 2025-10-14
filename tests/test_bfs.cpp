#include "catch_amalgamated.hpp"
#include "bfs.hpp"
#include "graph_utils.hpp"
#include <iostream>

using namespace std;

vector<pair<int, int>> edges_P = {{0,1}, {1,2}};
vector<pair<int, int>> edges_Q = {{0,1}, {1,2}, {2, 5}};
vector<pair<int, int>> edges_R = {{0,1}, {1,2}, {2, 5}, {4, 5}, {3, 4}, {0, 3}}; 

const int n = 3;

TEST_CASE("Shortest Path", "[3x3]") {
    auto adj_P = edge_list_to_adjacency_list(edges_P, n);
    auto adj_Q = edge_list_to_adjacency_list(edges_Q, n);
    auto adj_R = edge_list_to_adjacency_list(edges_R, n);
    auto g_P = path_to_edge_list(bfs(0, 2, adj_P, n));
    auto g_Q = path_to_edge_list(bfs(0, 2, adj_Q, n));
    auto g_R = path_to_edge_list(bfs(0, 2, adj_R, n));
    REQUIRE(g_P == edges_P);
    REQUIRE(g_Q != edges_Q);
    REQUIRE(g_Q == edges_P);
    REQUIRE(g_R != edges_R);
    //cout << vector_pair_to_string(g_P);
    //cout << vector_pair_to_string(g_Q);
    //cout << vector_pair_to_string(g_R);
    REQUIRE(g_R == edges_P);
}