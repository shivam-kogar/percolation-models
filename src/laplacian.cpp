#include <Eigen/Sparse>
#include <vector>
#include <utility>

using namespace std;
using namespace Eigen;

using SpMat = SparseMatrix<double>;
using TripletT = Triplet<double>;

/*
Construct Laplacian of a 2D n x n lattice with edges, using Eigen.
Assumes nodes are numbered from 0 to n^2 - 1.
*/
SpMat laplacian_2d(const vector< pair<int, int> >& edges, int n) {
    int N = n * n;
    vector<TripletT> triplets;
    vector<int> degree(N, 0);

    for (const auto& edge : edges) {
        int u = edge.first;
        int v = edge.second;

        triplets.emplace_back(u, v, -1.0);
        triplets.emplace_back(v, u, -1.0);
        degree[u]++;
        degree[v]++;
    }

    // Add diagonal entries (degrees)
    for (int i = 0; i < N; ++i) {
        triplets.emplace_back(i, i, static_cast<double>(degree[i]));
    }

    SpMat L(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}
