#include <iostream>
#include <iomanip>
#include <cmath>
#include "solver.hpp"
#include <fstream>
#include <filesystem>
#include <array>

struct LV_params {
    double alpha;
    double beta;
    double gamma;
    double delta;
};

// LABS EX1
Solver_ex1::Solver_ex1() {};


int Solver_ex1::func(double x, const double y[], double f[], void *params) {
    (void)(x);
    (void)(params);
    f[0] = (1) * y[0] * (1- y[0]);
    return GSL_SUCCESS;
}


void Solver_ex1::solve_ode(){
    std::filesystem::create_directories("data");
    // initialize_init_vals();

    const int n_loops = 1000;

    for (double x0 : init_values) {
        std::ofstream output("data/init_k1d_" + std::to_string(x0) + ".txt");
        
        output << "t" << "\t" << "x" <<  std::endl;
        std::cout << "=================\n";
        std::cout << "x(0) = " << x0 << std::endl;

        gsl_odeiv2_system sys = {func, nullptr, 1, nullptr};
        gsl_odeiv2_driver *driver =
            gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                          step, eps_abs, eps_rel);

        double t = t0;
        double y[1] = {x0};

        for (int j = 1; j <= n_loops; ++j) {
            double ti = t0 + j * (t1 - t0) / n_loops;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

            if (status != GSL_SUCCESS) {
                std::cerr << "error at x = " << ti << "\n";
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
// LABS EX2
Solver_ex2::Solver_ex2() {};

int Solver_ex2::func(double t, const double y[], double f[], void *params) {
    (void)(t); (void)(params);
    double x = y[0], yy = y[1];
    f[0] = x * (yy - 1.0);
    f[1] = 3.0 * x - 2.0 * yy + x*x - 2.0 * yy*yy;
    return GSL_SUCCESS;
}
struct Sys2Params { };
void Solver_ex2::write_vector_field_2d(double xmin, double xmax, int nx,
                                   double ymin, double ymax, int ny,
                                   const std::string &fname) {
    std::ofstream ofs(fname);
    ofs << "x\ty\tu\tv\n";
    for (int i = 0; i < nx; ++i) {
        double x = xmin + i * (xmax - xmin) / (nx - 1);
        for (int j = 0; j < ny; ++j) {
            double y = ymin + j * (ymax - ymin) / (ny - 1);
            double u = x * (y - 1.0);
            double v = 3.0 * x - 2.0 * y + x*x - 2.0 * y*y;
            ofs << x << "\t" << y << "\t" << u << "\t" << v << "\n";
        }
    }
    ofs.close();
}
void Solver_ex2::solve_ode(double t0, double t1, int n_steps, const std::string &out_prefix) {
    std::filesystem::create_directories("data");
    Sys2Params params;
    gsl_odeiv2_system sys = {func, nullptr, 2, &params};
    double abs_eps = 1e-8, rel_eps = 1e-8;
    double initial_step = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                      initial_step, abs_eps, rel_eps);

    std::vector<std::array<double,2> > initials = {{0.1, 0.1}, {0.5, 0.2}, {2.0, 1.5}, {-2.0, 0.5}, {-3.5, 1.0}, {0.0, -0.5}};

    for (auto ic : initials) {
        double x0 = ic[0], y0_val = ic[1];
        std::string fname = "data/" + out_prefix + "_x0_" + std::to_string(x0) + "_y0_" + std::to_string(y0_val) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\ty\n";

        double t = t0;
        double yvec[2] = {x0, y0_val};
        for (int j = 0; j <= n_steps; ++j) {
            double ti = t0 + j * (t1 - t0) / (double)n_steps;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, yvec);
            if (status != GSL_SUCCESS) { std::cerr << "GSL err\n"; break; }
            ofs << ti << "\t" << yvec[0] << "\t" << yvec[1] << "\n";
        }
        ofs.close();
    }
    gsl_odeiv2_driver_free(driver);
}


