#include "catch_amalgamated.hpp"
#include "laplacian.hpp"
#include "edge_list.hpp"
#include <armadillo>

using namespace arma;

auto edges_A = edge_list_2d(1, 3);
auto edges_B = edge_list_2d(0.5, 5);
auto edges_C = edge_list_2d(0.1, 10);

sp_mat A = laplacian_2d(edges_A, 3);
sp_mat B = laplacian_2d(edges_B, 5);
sp_mat C = laplacian_2d(edges_C, 10);

bool is_symmetric(sp_mat L) {
    return approx_equal(A, A.t(), "absdiff", 1e-12);
}

bool sum_to_zero(sp_mat L) {
    sp_vec row_sums = sum(L, 1); 
    sp_vec zeros(L.n_rows);                    

    return approx_equal(row_sums, zeros, "absdiff", 1e-12);
}

TEST_CASE("Symmetry of Laplacian") {
    REQUIRE(is_symmetric(A));
    REQUIRE(is_symmetric(B));
    REQUIRE(is_symmetric(C));
}

TEST_CASE("Rows and Columns Sum to Zero") {
    REQUIRE(sum_to_zero(A));
    REQUIRE(sum_to_zero(B));
    REQUIRE(sum_to_zero(C));
}


