# Duffing Oscillator Project

A modular C++ framework for simulating **classical and quantum dynamics of the Duffing oscillator**, with a focus on **quantum–classical transition studies** using the **Wigner function and Moyal equation**.

---

## Project Overview

This repository contains two main components:

1. **Classical Duffing oscillator analysis**
2. **Quantum phase-space dynamics via the Moyal equation**

The quantum solver evolves the Wigner distribution using a **Strang split-operator method** with FFT-based propagators.

---

## Directory Structure

### Classical Duffing Solver
include/
duffing.hpp
quantum_duffing.hpp
user_executable.hpp

src/
main.cpp
duffing.cpp
quantum_duffing.cpp

### Quantum Duffing Solver
moyal_solver/
include/
moyal_setup.hpp
moyal_solver.hpp
phase_space.hpp
potentials.hpp
propagators.hpp
wigner_df.hpp

src/
moyal_solver.cpp
phase_space.cpp
propagators.cpp
wigner_df.cpp

examples/
main.cpp

## Core Classes

- **`PhaseSpace`**  
  Manages phase-space grids, FFT plans, and coordinate transforms.

- **`WignerDistribution`**  
  Stores and manipulates the Wigner function, including initialization and diagnostics.

- **`Potential`**  
  Abstract base class for defining arbitrary time-dependent potentials.

- **`Propagator`**  
  Base class for kinetic and potential propagators used in operator splitting.

- **`MoyalSolver`**  
  Implements time evolution of the Wigner function using **Strang splitting**.

## Building

### Requirements
- C++17 compatible compiler
- CMake ≥ 3.12
- FFTW3

### Build Instructions
```bash
mkdir build
cd build
cmake ..
make