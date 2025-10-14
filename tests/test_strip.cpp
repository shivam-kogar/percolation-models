#include "catch_amalgamated.hpp"
#include "strip_utils.hpp"
#include "graph_utils.hpp"
#include "bfs.hpp"
#include <iostream>
#include <algorithm>

using namespace std;

bool strip_correct(const vector<pair<int, int>>& strip, pair<int, int> edge);

const int n = 3;

vector<pair<int, int>> edge_S = {{0,1}, {1,2}, {0,3}, {1,4}, {2,5}, {3,4}, {4,5}, {6,7}, {7,8}, {3,6}, {4, 7}, {5, 8}};

auto adj_S = edge_list_to_adjacency_list(edge_S, n); 

auto geodesic = bfs(0, 2, adj_S, n);
auto geodesic_2 = bfs(3, 5, adj_S, n);

auto strip_1 = get_strip_edges_from_ids(geodesic, edge_S, 1, n, true);
auto strip_2 = get_strip_edges_from_ids(geodesic, edge_S, 2, n, true);
auto strip_3 = get_strip_edges_from_ids(geodesic_2, edge_S, 1, n, true);
auto strip_4 = get_strip_edges_from_ids(geodesic_2, edge_S, 0, n, true);


TEST_CASE("Strip Edges") {
    cout << vector_pair_to_string(strip_1);
    REQUIRE((strip_1, {0,1}));
    REQUIRE(!strip_has_edge(strip_1, {0,2}));
    REQUIRE(strip_has_edge(strip_1, {1,2}));
    REQUIRE(strip_has_edge(strip_1, {0,3}));
    REQUIRE(!strip_has_edge(strip_1, {5,8}));
    REQUIRE(strip_has_edge(strip_2, {5,8}));
    REQUIRE(strip_has_edge(strip_1, {3,4}));
    REQUIRE(strip_has_edge(strip_1, {4,5}));
    REQUIRE(strip_has_edge(strip_1, {2,5}));
    REQUIRE(strip_has_edge(strip_1, {1,4}));
    REQUIRE(!strip_has_edge(strip_1, {6,7}));
    REQUIRE(strip_has_edge(strip_2, {6,7}));
    REQUIRE(strip_has_edge(strip_2, {4,7}));
    for (auto edge : edge_S) {
        REQUIRE(strip_has_edge(strip_3, edge));
    }
    auto g2_edges = path_to_edge_list(geodesic_2);
    for (auto edge : edge_S) {
        if (strip_has_edge(g2_edges, edge)) {
            REQUIRE(strip_has_edge(strip_4, edge));
        }
        else {
            REQUIRE(!strip_has_edge(strip_4, edge));
        }
    }
}

bool strip_has_edge(const vector<pair<int, int>>& strip, pair<int, int> edge) {
    return find(strip.begin(), strip.end(), edge) != strip.end();
}