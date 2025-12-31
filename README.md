# Transport in 2D Bond Percolation Networks

This repository contains C++ simulations used to measure transport in
2D bond percolation networks. We generate a square lattice parameterized
by a bond probability 'p', which is the probability of a bond existing
between any two lattice nodes. We fix endpoints distant from the lattice
boundary and, if some path exists connecting them, we measure the
- total resistance between the endpoints (using Kirchhoff's and Ohm's laws
  and considering the whole network)
- resistance along the geodesic (i.e., the length of the geodesic path between
  the endpoints)
- resistance within a strip about the geodesic
And we take the reciprocal of each quantity to obtain network conductances,
geodesic conductances, and strip conductances, respectively.


