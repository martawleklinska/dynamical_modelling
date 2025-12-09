#include "bifurcation.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <array>

Solver_G::Solver_G(double m_) : m(m_) {}

int Solver_G::func(double t, const double y[], double f[], void *params) {
    (void)(t); (void)(params);
    double x = y[0];
    double m = 0.0;
    f[0] = m - std::abs(x);
    return GSL_SUCCESS;
}

void Solver_G::solve() {
    std::filesystem::create_directories("data");
    gsl_odeiv2_system sys = {func, nullptr, 1, nullptr};

    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, eps_abs, eps_rel);

    std::vector<double> initials;
    for (double x0 = -2.0; x0 <= 2.0; x0 += 0.2) initials.push_back(x0);

    for (auto x0 : initials) {
        std::string fname = "data/ex1_m" + std::to_string(m) + "_x0_" + std::to_string(x0) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\n";

        double t = t0;
        double y[1] = {x0};

        for (int i = 0; i <= n_steps; ++i) {
            double ti = t0 + i * (t1 - t0) / n_steps;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
            if (status != GSL_SUCCESS) break;
            ofs << ti << "\t" << y[0] << "\n";
        }
        ofs.close();
    }
    gsl_odeiv2_driver_free(driver);
}

// ===================== H system =====================

Solver_H::Solver_H(double p_) : p(p_) {}

int Solver_H::func(double t, const double y[], double f[], void *params) {
    (void)(t); (void)(params);
    double x = y[0];
    double p = 5.;
    f[0] = (x - 1) * (x * x + 2 * x - p);
    return GSL_SUCCESS;
}

void Solver_H::solve() {
    std::filesystem::create_directories("data");
    gsl_odeiv2_system sys = {func, nullptr, 1, nullptr};

    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, eps_abs, eps_rel);

    std::vector<double> initials;
    for (double x0 = -4.0; x0 <= 2.0; x0 += 0.2) initials.push_back(x0);

    for (auto x0 : initials) {
        std::string fname = "data/ex2_" + std::to_string(p) + "_x0_" + std::to_string(x0) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\n";

        double t = t0;
        double y[1] = {x0};

        for (int i = 0; i <= n_steps; ++i) {
            double ti = t0 + i * (t1 - t0) / n_steps;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
            if (status != GSL_SUCCESS) break;
            ofs << ti << "\t" << y[0] << "\n";
        }
        ofs.close();
    }
    gsl_odeiv2_driver_free(driver);
}
