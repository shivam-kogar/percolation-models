#pragma once
#include <Eigen/Sparse>
#include <vector>
#include <utility>

using namespace std;
using namespace Eigen;

using SpMat = Eigen::SparseMatrix<double>;

SpMat laplacian_2d(const vector<pair<int, int>>& edges, int n);
