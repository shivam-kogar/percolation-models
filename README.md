# Transport in 2D Bond Percolation Networks

This repository contains a **C++ implementation for simulating transport in 2D bond percolation networks**.  
The program generates square lattices with a given bond probability and computes transport-related quantities on the resulting random network.

The focus of this codebase is **numerical simulation and data generation**. Results are written to **CSV files** for external analysis. No visualization is performed by the code itself.

---

## Features

- Construction of 2D square-lattice bond percolation networks
- Configurable bond probability and lattice parameters (defined in code)
- Computation of transport-related observables
- CSV output suitable for post-processing
- Modular C++ design using separate source and header files
- Uses the Eigen library for linear algebra

---

## Repository Structure
- include/ # Header files
- src/ # Source code modules
- tests/ # Test code
- latex-notes/ # Notes and documentation (LaTeX)
- main.cpp # Program entry point
- README.md # This file

---

## Requirements

- A C++ compiler supporting **C++17 or later**
- Tested with `g++-15`
- Eigen (header-only linear algebra library)

---

## Build Instructions

This project is compiled using `g++` with C++17 and Eigen.

### macOS (Homebrew GCC)

```bash
g++-15 -std=c++17 \
  src/*.cpp main.cpp \
  -Iinclude \
  -I/opt/homebrew/include/eigen3 \
  -O2 \
  -o percolation
