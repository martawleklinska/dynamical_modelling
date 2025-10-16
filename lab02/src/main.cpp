#include <iostream>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "solver.hpp"

int main() {
    Solver sol;

    std::cout << "test" << std::endl;

    // sol.solve_with_rk8pd();
    // std::cout << "done pd8" << std::endl;

    // ex1: solving the \dot{x}=kx(1-x)
    // Solver_ex1 sol_ex1;
    // sol_ex1.solve_ode();

    Solver_ex2 sol_ex2;
    sol_ex2.write_vector_field_2d(-6.0, 6.0, 31, -3.0, 4.0, 25, "data/sys2_vector_field.txt");


    sol_ex2.solve_ode(0.0, 100.0, 5000, "sys2_traj");
    std::cout << "saved ex2 data/.\n";


    return 0;
}
