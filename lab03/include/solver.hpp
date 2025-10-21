#pragma once
#include <vector>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


class Solver_pendulum {
public:
    Solver_pendulum(int model_type); 
    void solve(double t0, double t1, int n_steps);

    void analyze_fixed_points();

private:
    static int func(double t, const double y[], double f[], void *params);
    static void jacobian(double theta, double omega, gsl_matrix *J, int model_type);

    int model; // 1,2,3 oznacza równania (3),(4),(5)
    double abs_eps = 1e-8;
    double rel_eps = 1e-8;
};