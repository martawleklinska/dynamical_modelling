#pragma once
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

class Solver {
public:
    Solver();
    void solve_with_rk4();
    void solve_with_rk8pd();

private:
    static int func(double x, const double y[], double f[], void *params);
    static int jac(double x, const double y[], double *dfdy, double dfdt[], void *params);
    double y_exact(double x, double y0);

    double x0 = 0.0;
    double x1 = 1.0;
    double step = 1e-3;
    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    std::vector<double> init_values = {0.0, 0.5, 1.0, 2.0};
};
