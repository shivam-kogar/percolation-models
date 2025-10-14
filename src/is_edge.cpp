#include <iostream>
#include <random>

using namespace std;

/*
Uses percolation constant p to 'decide' whether to construct an edge between two points on the lattice via a 
pseudo-random number generator (PRNG)
*/
bool is_edge(double p) {
    mt19937 gen(random_device{}()); // PRNG from std library
    uniform_real_distribution<double> dist(0.0, 1.0); // generates random number between 0 and 1

    double randNum = dist(gen);
    return randNum < p; // return true if the random number is less than p, return false else
}