#pragma once
#include <Eigen/Sparse>

using SpMat = Eigen::SparseMatrix<double>;

double resistance_distance_2d(int i, int j, const SpMat& L, double s);
