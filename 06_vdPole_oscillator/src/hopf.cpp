#include "hopf.hpp"
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iostream>

int SolverVdP::func(double t, const double y[], double f[], void *params){
    (void)(t);
    auto *pp = static_cast<VdPParams*>(params);
    double eps = pp->eps;
    double m = pp->m;
    double x = y[0];
    double u = y[1];
    f[0] = u;
    f[1] = -eps*(x*x - 1.0)*u - x + m;
    return GSL_SUCCESS;
}

void SolverVdP::solve_trajectories(const std::vector<std::array<double,2>> &initials,
                                   const std::string &out_prefix)
{
    std::filesystem::create_directories("data");
    gsl_odeiv2_system sys = {func, nullptr, 2, &p};
    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-9, 1e-9);

    for (auto ic : initials) {
        double x0 = ic[0], u0 = ic[1];
        std::string fname = "data/" + out_prefix + "_x0_" + std::to_string(x0) + "_u0_" + std::to_string(u0) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\tu\n";
        double t = t0;
        double y[2] = {x0, u0};
        for (int i = 0; i <= n_steps; ++i) {
            double ti = t0 + i * (t1 - t0) / double(n_steps);
            int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
            if (status != GSL_SUCCESS) { std::cerr << "GSL error\n"; break; }
            ofs << ti << "\t" << y[0] << "\t" << y[1] << "\n";
        }
        ofs.close();
    }
    gsl_odeiv2_driver_free(d);
}
