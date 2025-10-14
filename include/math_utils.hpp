#pragma once
#include <vector>
#include <numeric>
using namespace std;

inline double mean(const vector<double>& data) {
    if (data.empty()) return 0.0;
    return accumulate(data.begin(), data.end(), 0.0) / data.size();
}

inline double stdev(const vector<double>& data, double mean) {
    double sum_sq = 0.0;
    for (double x : data) {
        sum_sq += (x - mean) * (x - mean);
    }
    return sqrt(sum_sq / (data.size() - 1));
}