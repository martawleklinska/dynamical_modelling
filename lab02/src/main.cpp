#include <iostream>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "solver.hpp"

int main() {
    Solver sol;

    std::cout << "test" << std::endl;

    sol.solve_with_rk4();
    std::cout << "done" << std::endl;

    sol.solve_with_rk8pd();
    std::cout << "done pd8" << std::endl;

    return 0;
}
