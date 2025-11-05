#pragma once
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

struct DuffingParams {
    double zeta, alpha, beta, gamma, omega;
};

class Duffing {
private:
    static int func(double t, const double y[], double f[], void *params);
    DuffingParams p; 

    double t0 = 0.0;
    double t1 = 10.;
    int n_steps = 10000;

public:
    Duffing(double zeta = .05, double alpha = 1., double beta = 1., 
            double gamma = .2, double omega = 1.);
    void solve();
};
