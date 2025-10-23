#include "solver.hpp"
#include <iostream>

int main() {

    for (int model = 1; model <= 3; ++model) {
        std::cout << "\n=== Model " << model << " ===\n";
        Solver_pendulum solver(model);
        solver.solve(0.0, 10.0, 20000);
    }
    std::cout << "Wyniki zapisane w folderze data/.\n";
    return 0;
}
