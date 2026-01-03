# Duffing Oscillator Project

A C++ project for simulating **classical and quantum dynamics of the Duffing oscillator**, including
nonlinear dynamics analysis and **quantum phase-space evolution using the Wigner function and Moyal equation**.

The project consists of **two independent executables**:
1. A **classical Duffing solver** with an interactive analysis menu
2. A **quantum Wigner/Moyal solver** for phase-space dynamics

---

## Classical Duffing Solver

### Build and Run

From the project root:

```bash
cd projekt_duffing/
mkdir build
cd build
cmake ..
make
./main
```
This launches an interactive menu:
```bash
====== DUFFING OSCILLATOR SOLVER ======

=========== DUFFING OSCILLATOR MENU ===========
1. ODE solver (multiple parameter sets)
2. Poincaré map analysis
3. Bifurcation analysis
4. Lyapunov exponent analysis
5. Energy analysis
6. Quantum Duffing analysis
0. Exit
===============================================
Please enter your choice (0-6):
```
### Changing System Parameters

All classical Duffing parameters (drive strength, damping, nonlinearity, frequency, etc.)
are defined in:
```bash
src/user_executable.cpp
```
## Quantum Wigner / Moyal Solver

This part of the project evolves the Wigner distribution function using the Moyal equation and a Strang split-operator method with FFTs.

### Build and Run
```bash
cd projekt_duffing/moyal_solver/
./build.sh 
cd build  
./main
```
It initializes a Gaussian (coherent) state in phase space;

Evolves it under a chosen potential (e.g. Duffing, but you can chose another in main.cpp);

Outputs Wigner function data files to build/output/;

Visualization done in Julia (main.jl);

Simulation parameters (grid size, time step, effective ħ, initial state) are set in the example/main.cpp.

## Requirements

C++17 compatible compiler

CMake ≥ 3.12

FFTW3 (required for quantum solver)

GSL