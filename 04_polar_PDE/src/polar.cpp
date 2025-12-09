#include "polar.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>

Polar_solver::Polar_solver() {}

int Polar_solver::func(double t, const double y[], double f[], void *params) {
    (void)(t); (void)(params);
    double x = y[0];
    double yy = y[1];

    double r2 = x*x + yy*yy;
    double F = (r2*r2 - 4.0*r2 + 1.0); 

    f[0] = -2.0 * yy - x * F;
    f[1] =  2.0 * x  - yy * F;
    return GSL_SUCCESS;
}

void Polar_solver::solve_ode(double t0, double t1, int n_steps, const std::string &out_prefix) {
    std::filesystem::create_directories("data");

    gsl_odeiv2_system sys = {func, nullptr, 2, nullptr};

    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    double h = (t1 - t0) / n_steps;

    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, eps_abs, eps_rel);

    std::vector<double> x0_list = {0.2, 0.4, 0.5175, 0.5180, 1.0, 1.8, 2.5};
    double y0 = 0.0;

    for (double x0 : x0_list) {
        std::string fname = "data/" + out_prefix + "_x0_" + std::to_string(x0) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\ty\n";

        double yvec[2] = {x0, y0};
        double t = t0;

        for (int i = 0; i <= n_steps; ++i) {
            double ti = t0 + i * (t1 - t0) / n_steps;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, yvec);
            if (status != GSL_SUCCESS) {
                std::cerr << "GSL error at step " << i << "\n";
                break;
            }
            ofs << ti << "\t" << yvec[0] << "\t" << yvec[1] << "\n";
        }

        ofs.close();
    }

    gsl_odeiv2_driver_free(driver);
}


// ====================== radius ===============
int Polar_solver::func_rad(double x, const double y[], double f[], void *params) {
    (void)(x);
    (void)(params);
    f[0] = -y[0] * (pow(y[0], 4)-4 * y[0] * y[0] + 1);
    return GSL_SUCCESS;
}


void Polar_solver::solve_radius_ode(double t0, double t1, int n_steps, const std::string &out_prefix){
    std::filesystem::create_directories("data");
    std::vector<double> r0_list = {0.2, 0.4, 0.5175, 0.5180, 1.0, 1.8, 2.5};

    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    double h = (t1 - t0) / n_steps;

    for (double x0 : r0_list) {
        std::ofstream output("data/" + out_prefix + "_x0_" + std::to_string(x0) + ".txt");
        
        output << "t" << "\t" << "r" <<  std::endl;
        std::cout << "=================\n";
        std::cout << "r(0) = " << x0 << std::endl;

        gsl_odeiv2_system sys = {func_rad, nullptr, 1, nullptr};
        gsl_odeiv2_driver *driver =
            gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                          h, eps_abs, eps_rel);

        double t = t0;
        double y[1] = {x0};

        for (int j = 1; j <= n_steps; ++j) {
            double ti = t0 + j * (t1 - t0) / n_steps;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

            if (status != GSL_SUCCESS) {
                std::cerr << "error at r = " << ti << "\n";
                break;
            }

            output << ti << "\t" << y[0] << "\n";

            if (std::abs(y[0]) > 1e6) {
                std::cout << "t " << ti << "\n";
                break;
            }
        }

        gsl_odeiv2_driver_free(driver);
        output.close();
    }

}
