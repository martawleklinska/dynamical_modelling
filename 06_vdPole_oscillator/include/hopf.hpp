#pragma once
#include <vector>
#include <array>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

struct VdPParams { double eps; double m; };

class SolverVdP {
private:
    static int func(double t, const double y[], double f[], void *params);
    VdPParams p;
    double t0 = 0.0;
    double t1 = 200.0;  
    int n_steps = 20000;
public:
    SolverVdP(double eps=0.5, double m=0.0) { p.eps = eps; p.m = m; }
    void set_params(double eps_, double m_) { p.eps = eps_; p.m = m_; }
    void solve_trajectories(const std::vector<std::array<double,2>> &initials,
                            const std::string &out_prefix);
};
