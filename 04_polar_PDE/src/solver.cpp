#include"solver.hpp"
#include<iostream>
#include<fstream>
#include<cmath>
#include<array>
#include<filesystem>

Solver::Solver(){};

int Solver::func(double t, const double y[], double f[], void *params){
    (void)(t); (void)(params);
    double x = y[0];
    double u = y[1];
    f[0] = u;
    f[1] = - (2.0/3.0) * u * u * u + 2.0 * u - x;
    return GSL_SUCCESS;
}

void Solver::solve_ode(){
    std::filesystem::create_directories("data");

    gsl_odeiv2_system sys = {func, nullptr, 2, nullptr};

    double eps_abs = 1e-08;
    double eps_rel = 1e-08;
    double h = (t1-t0)/n_steps;

    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, eps_abs, eps_rel);
    
    std::vector<std::array<double,2> > initials = {{0.0, 4.0}, {0.0, -4.0}, {-0.7, -0.5}, {-1.5, -3}, {1., -3.}, {0.5, -3.}};
    for (auto ic : initials){
        std::string fname = "data/ex2_x0_" + std::to_string(ic[0]) + "_u0_" + std::to_string(ic[1])+ ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\tu\n";
        double x0 = ic[0], u0_val = ic[1];
        double t = t0;
        double yvec[2] = {x0, u0_val};
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

