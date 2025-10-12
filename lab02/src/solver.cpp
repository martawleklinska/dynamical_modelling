#include <iostream>
#include <iomanip>
#include <cmath>
#include "solver.hpp"
#include <fstream>
#include <filesystem>

// structure to hold parameters
struct LV_params {
    double alpha;
    double beta;
    double gamma;
    double delta;
};

Solver::Solver() {}

// \dot{x} = \alpha * x - \beta * x * y
// \dot{y} = - \gamma * y + \delta * x * y
int Solver::func(double t, const double y[], double f[], void *params) {
    (void)(t);
    LV_params *p = static_cast<LV_params *>(params);
    double alpha = p->alpha;
    double beta  = p->beta;
    double gamma = p->gamma;
    double delta = p->delta;

    f[0] = y[0] * (alpha - beta * y[1]);
    f[1] = y[1] * (-gamma + delta * y[0]);

    return GSL_SUCCESS;
}

// Jacobian matrix:
// df/dy = [ α - βy,   -βx
//            δy,     -γ + δx ]
int Solver::jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    (void)(t);
    LV_params *p = static_cast<LV_params *>(params);
    double alpha = p->alpha;
    double beta  = p->beta;
    double gamma = p->gamma;
    double delta = p->delta;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, alpha - beta * y[1]);
    gsl_matrix_set(m, 0, 1, -beta * y[0]);
    gsl_matrix_set(m, 1, 0, delta * y[1]);
    gsl_matrix_set(m, 1, 1, -gamma + delta * y[0]);

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;

    return GSL_SUCCESS;
}

void Solver::solve_with_rk4() {
    std::filesystem::create_directories("data");

    // define parameters
    LV_params params = {1.0, 0.5, 2.0, 0.3};

    const int n_loops = 1000;

    for (double y0 : init_values) {
        std::ofstream output("data/lv_rk4_y0_" + std::to_string(y0) + ".txt");
        output << "t\tx\ty\n";
        std::cout << "=================\n";
        std::cout << "Initial x(0)=1.0, y(0)=" << y0 << std::endl;

        gsl_odeiv2_system sys = {func, jac, 2, &params};
        gsl_odeiv2_driver *driver =
            gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4,
                                          step, eps_abs, eps_rel);

        double t = x0;
        double y[2] = {1.0, y0}; // initial populations

        for (int j = 1; j <= n_loops; ++j) {
            double ti = x0 + j * (x1 - x0) / n_loops;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

            if (status != GSL_SUCCESS) {
                std::cerr << "error at t = " << ti << "\n";
                break;
            }

            output << ti << "\t" << y[0] << "\t" << y[1] << "\n";

            if (std::abs(y[0]) > 1e6 || std::abs(y[1]) > 1e6) {
                std::cout << "Blow-up near t ≈ " << ti << "\n";
                break;
            }
        }

        gsl_odeiv2_driver_free(driver);
        output.close();
    }
}

void Solver::solve_with_rk8pd() {
    std::filesystem::create_directories("data");

    LV_params params = {1.0, 0.5, 2.0, 0.3};
    const int n_loops = 1000;

    for (double y0 : init_values) {
        std::ofstream output("data/lv_rk8pd_y0_" + std::to_string(y0) + ".txt");
        output << "t\tx\ty\n";
        std::cout << "=================\n";
        std::cout << "Initial x(0)=1.0, y(0)=" << y0 << std::endl;

        gsl_odeiv2_system sys = {func, jac, 2, &params};
        gsl_odeiv2_driver *driver =
            gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                          step, eps_abs, eps_rel);

        double t = x0;
        double y[2] = {1.0, y0};

        for (int j = 1; j <= n_loops; ++j) {
            double ti = x0 + j * (x1 - x0) / n_loops;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

            if (status != GSL_SUCCESS) {
                std::cerr << "error at t = " << ti << "\n";
                break;
            }

            output << ti << "\t" << y[0] << "\t" << y[1] << "\n";

            if (std::abs(y[0]) > 1e6 || std::abs(y[1]) > 1e6) {
                std::cout << "Blow-up near t ≈ " << ti << "\n";
                break;
            }
        }

        gsl_odeiv2_driver_free(driver);
        output.close();
    }
}
