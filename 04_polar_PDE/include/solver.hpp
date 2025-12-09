#pragma once
#include <vector>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

class Solver{
    private:
        static int func(double t, const double y[], double f[], void *params);
        double t0 = 0.0;
        double t1 = 60.;
        int n_steps = 5000;
    public:
        Solver();
        void solve_ode();
};