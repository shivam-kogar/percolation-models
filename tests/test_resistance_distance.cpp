#include "catch_amalgamated.hpp"
#include "laplacian.hpp"
#include "resistance_distance.hpp"
#include <armadillo>
#include <iostream>

using namespace std;

vector<pair<int, int>> edges_X = {{0,1}, {1,2}};
vector<pair<int, int>> edges_Y = {{0,1}, {1,2}, {2, 5}};
vector<pair<int, int>> edges_Z = {{0,1}, {1,2}, {2, 5}, {4, 5}, {3, 4}, {0, 3}}; 

sp_mat X = laplacian_2d(edges_X, 3);
sp_mat Y = laplacian_2d(edges_Y, 3);
sp_mat Z = laplacian_2d(edges_Z, 3);

bool rd_correct(int i, int j, sp_mat L, double guess) {
    // cout << "\n" << resistance_distance_2d(i, j, L, 1e-12) << "\n";
    return abs(resistance_distance_2d(i, j, L, 1e-12) - guess) < 1e-4;
}

TEST_CASE("Resistance Distance") {
    REQUIRE(rd_correct(0, 2, X, 2));
    REQUIRE(rd_correct(0, 2, Y, 2));
    REQUIRE(!rd_correct(0, 2, Z, 2));
    REQUIRE(rd_correct(0, 2, Z, 1.3333));
}




