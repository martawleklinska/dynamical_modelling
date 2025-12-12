#pragma once
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <string>
#include<vector>
#include<array>

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
    std::vector<std::array<double, 2>> initials = {
        {0.5, -2.0}, {-1.0, 2.0}, {-1., 0.50}, {0.5, -1.70},
        {0.0, 0.20}, {0.0, 0.50}};

public:
    Duffing(double zeta = .05, double alpha = 1., double beta = 1., 
            double gamma = .2, double omega = 1.);
    void solve(std::string filename);

    void poincare_map(double discard_transient,
                      int n_periods_sample,
                      const std::string &out_prefix,
                      double t0_override = -1.0);

    void bifurcation_scan(double gamma_min, double gamma_max, int n_gamma,
                          double discard_transient, int samples_per_gamma,
                          const std::string &out_prefix,
                          double t0_override = -1.0);

};
