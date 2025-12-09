#pragma once
#include <vector>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

class Solver_G {
private:
    static int func(double t, const double y[], double f[], void *params);
    double m = 0.0;
    double t0 = 0.0;
    double t1 = 10.0;
    int n_steps = 2000;

public:
    Solver_G(double m_);
    void solve();
};

class Solver_H {
private:
    static int func(double t, const double y[], double f[], void *params);
    double p = 5.;
    double t0 = 0.0;
    double t1 = 2.0;
    int n_steps = 2000;

public:
    Solver_H(double p_);
    void solve();
};
