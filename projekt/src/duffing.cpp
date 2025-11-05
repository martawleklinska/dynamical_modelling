#include "duffing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <array>

Duffing::Duffing(double zeta_, double alpha_, double beta_, double gamma_, double omega_) {
    p = {zeta_, alpha_, beta_, gamma_, omega_};
}

int Duffing::func(double t, const double y[], double f[], void *params) {
    auto *p = static_cast<DuffingParams*>(params);
    double x = y[0];
    double v = y[1];

    f[0] = v;
    f[1] = -2.0 * p->zeta * v + p->alpha * x - p->beta * std::pow(x, 3) + p->gamma * std::cos(p->omega * t);
    return GSL_SUCCESS;
}

void Duffing::solve() {
    std::filesystem::create_directories("data");

    gsl_odeiv2_system sys = {func, nullptr, 2, &p};

    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, h, 1e-8, 1e-8);

    std::vector<std::array<double, 2>> initials = {
        {0.5, -2.}, {-1.0, 2.0}, {-1., 0.5}, {0.5, -1.7}};

    for (auto ic : initials) {
        std::string fname = "data/duff_x0_" + std::to_string(ic[0]) +
                            "_v0_" + std::to_string(ic[1]) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\tv\n";

        double y[2] = {ic[0], ic[1]};
        double t = t0;

        for (int i = 0; i <= n_steps; ++i) {
            double ti = t0 + i * (t1 - t0) / n_steps;
            int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
            if (status != GSL_SUCCESS) {
                std::cerr << "GSL error at step " << i << "\n";
                break;
            }
            ofs << ti << "\t" << y[0] << "\t" << y[1] << "\n";
        }
        ofs.close();
    }

    gsl_odeiv2_driver_free(d);
}
