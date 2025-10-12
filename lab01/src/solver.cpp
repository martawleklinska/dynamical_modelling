#include <iostream>
#include <iomanip>
#include <cmath>
#include "solver.hpp"
#include <fstream>
#include<filesystem>

Solver::Solver() {}

int Solver::func(double x, const double y[], double f[], void *params) {
    (void)(x);
    (void)(params);
    f[0] = y[0] * y[0];
    return GSL_SUCCESS;
}

double Solver::y_exact(double x, double y0){
    return 1.0 / ((1.0/y0) - x);
}

int Solver::jac(double x, const double y[], double *dfdy, double dfdt[], void *params) {
    (void)(x);
    (void)(params);
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
    gsl_matrix_set(&dfdy_mat.matrix, 0, 0, 2.0 * y[0]);
    dfdt[0] = 0.0;
    return GSL_SUCCESS;
}

void Solver::solve_with_rk4() {
    std::filesystem::create_directories("data");

    const int n_loops = 1000;

    for (double y0 : init_values) {
        std::ofstream output("data/init_cond_" + std::to_string(y0) + ".txt");
        std::ofstream output_err("data/error_" + std::to_string(y0) + ".txt");
        std::ofstream output_errrel("data/errorrel_" + std::to_string(y0) + ".txt");
        
        output << "x" << "\t" << "y" <<  std::endl;
        std::cout << "=================\n";
        std::cout << "y(0) = " << y0 << std::endl;

        gsl_odeiv2_system sys = {func, jac, 1, nullptr};
        gsl_odeiv2_driver *driver =
            gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4,
                                          step, eps_abs, eps_rel);

        double x = x0;
        double y[1] = {y0};

        for (int j = 1; j <= n_loops; ++j) {
            double xi = x0 + j * (x1 - x0) / n_loops;
            int status = gsl_odeiv2_driver_apply(driver, &x, xi, y);

            if (status != GSL_SUCCESS) {
                std::cerr << "error at x = " << xi << "\n";
                break;
            }

            double yexa = y_exact(xi, y0);

            output << xi << "\t" << y[0] << "\n";
            output_err << xi << "\t" << std::abs(y[0] - yexa) << "\n";
            output_errrel << xi << "\t" << std::abs(y[0] - yexa)/yexa << "\n";

            if (std::abs(y[0]) > 1e6) {
                std::cout << "x " << xi << "\n";
                break;
            }
        }

        gsl_odeiv2_driver_free(driver);
        output.close();
        output_err.close();
        output_errrel.close();
    }
}


void Solver::solve_with_rk8pd() {
    std::filesystem::create_directories("data");

    const int n_loops = 1000;

    for (double y0 : init_values) {
        std::ofstream output("data/init_cond8_" + std::to_string(y0) + ".txt");
        std::ofstream output_err("data/error8_" + std::to_string(y0) + ".txt");
        std::ofstream output_errrel("data/errorrel8_" + std::to_string(y0) + ".txt");
        
        output << "x" << "\t" << "y" <<  std::endl;
        std::cout << "=================\n";
        std::cout << "y(0) = " << y0 << std::endl;

        gsl_odeiv2_system sys = {func, jac, 1, nullptr};
        gsl_odeiv2_driver *driver =
            gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                          step, eps_abs, eps_rel);

        double x = x0;
        double y[1] = {y0};

        for (int j = 1; j <= n_loops; ++j) {
            double xi = x0 + j * (x1 - x0) / n_loops;
            int status = gsl_odeiv2_driver_apply(driver, &x, xi, y);

            if (status != GSL_SUCCESS) {
                std::cerr << "error at x = " << xi << "\n";
                break;
            }

            double yexa = y_exact(xi, y0);

            output << xi << "\t" << y[0] << "\n";
            output_err << xi << "\t" << std::abs(y[0] - yexa) << "\n";
            output_errrel << xi << "\t" << std::abs(y[0] - yexa)/yexa << "\n";

            if (std::abs(y[0]) > 1e6) {
                std::cout << "x " << xi << "\n";
                break;
            }
        }

        gsl_odeiv2_driver_free(driver);
        output.close();
        output_err.close();
        output_errrel.close();
    }
}