#define CATCH_CONFIG_MAIN
#include "catch_amalgamated.hpp"
#include "edge_list.hpp"

TEST_CASE("Size of Edge List") {
    REQUIRE(edge_list_2d(1, 2).size() == 4);
    REQUIRE(edge_list_2d(1, 3).size() == 12);
    REQUIRE(edge_list_2d(0, 5).size() == 0);
    REQUIRE(edge_list_2d(0, 10).size() == 0);
}
