#pragma once
#include <vector>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


class Polar_solver {
public:
    Polar_solver();
    void solve_ode(double t0, double t1, int n_steps, const std::string &out_prefix);
    void solve_radius_ode(double t0, double t1, int n_steps, const std::string &out_prefix);

private:
    static int func(double t, const double y[], double f[], void *params);
    static int func_rad(double x, const double y[], double f[], void *params);
};