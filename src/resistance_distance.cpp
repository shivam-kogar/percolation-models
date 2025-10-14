#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <limits>

using SpMat = Eigen::SparseMatrix<double>;
using Vec   = Eigen::VectorXd;

double resistance_distance_2d(int i, int j, const SpMat& L, double s) {
    SpMat A = L;
    const int N = A.rows();
    for (int k = 0; k < N; ++k) {
        A.coeffRef(k, k) += s;
    }

    // Set up Conjugate Gradient solver with diagonal preconditioner
    Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> solver;
    solver.setTolerance(1e-8);
    solver.setMaxIterations(500);
    A.makeCompressed();
    solver.compute(A);

    if (solver.info() != Eigen::Success) {
        return std::numeric_limits<double>::infinity();
    }

    Vec e_j = Vec::Zero(N);
    e_j(j) = 1.0;
    Vec x_j = solver.solve(e_j);
    if (solver.info() != Eigen::Success) {
        return std::numeric_limits<double>::infinity();
    }

    Vec e_i = Vec::Zero(N);
    e_i(i) = 1.0;
    Vec x_i = solver.solve(e_i);
    if (solver.info() != Eigen::Success) {
        return std::numeric_limits<double>::infinity();
    }

    double R_ij = x_i(i) + x_j(j) - x_i(j) - x_j(i);
    return R_ij;
}
