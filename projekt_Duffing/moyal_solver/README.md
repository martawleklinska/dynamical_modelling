# Structured Moyal Equation Solver

A modern, object-oriented C++ implementation of the Moyal equation solver using Strang splitting method for quantum phase space dynamics.

## Features

- **Object-oriented design**: Separate classes for phase space, Wigner distribution, propagators, and solver
- **Flexible potential system**: Easy to define custom potentials
- **FFTW integration**: Efficient FFT operations for split-operator method
- **Configurable parameters**: No more hardcoded values
- **Clean API**: Easy to use and extend

## Class Structure

- `PhaseSpace`: Manages coordinate grids and FFT operations
- `WignerDistribution`: Handles Wigner function operations and analysis
- `Propagator`: Base class for kinetic and potential propagators
- `MoyalSolver`: Main solver class implementing Strang splitting
- `Potential`: Abstract base for defining potential functions

## Building

Requirements:
- C++17 compatible compiler
- CMake 3.12+
- FFTW3 library

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

```cpp
#include "moyal_solver.hpp"
#include "potentials.hpp"

// Configure the solver
MoyalConfig config;
config.gridX = 1024;
config.gridP = 512;
config.timeStepsN = 800;
// ... set other parameters

// Create potential
auto potential = std::make_unique<GaussianPotential>(amplitude, center, width);

// Create and run solver
MoyalSolver solver(config, std::move(potential));
solver.setInitialGaussian(x0, p0, sigma_x, sigma_p);
solver.evolve();
```

## Comparison with Original Code

### Original Structure
- Single monolithic file (`s.cpp`, ~1600 lines)
- Global variables and hardcoded parameters
- Manual memory management
- Mixed concerns (FFT, physics, I/O all together)

### Structured Implementation
- Modular design with clear separation of concerns
- Configuration-driven parameters
- Automatic memory management with RAII
- Extensible potential system
- Clean, documented API

## Physics Background

This solver implements the Moyal equation (quantum Liouville equation) in phase space using the Wigner distribution function. The time evolution uses Strang splitting:

```
exp(-i*dt*H) ≈ exp(-i*dt/2*T) exp(-i*dt*V) exp(-i*dt/2*T)
```

where T is kinetic energy and V is potential energy operator.

## Examples

See `examples/basic_example.cpp` for a complete working example that reproduces the behavior of the original solver with the structured approach.
