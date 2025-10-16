#pragma once
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

class Solver {
public:
    Solver();
    void solve_with_rk8pd();

private:
    static int func(double x, const double y[], double f[], void *params);
    static int jac(double x, const double y[], double *dfdy, double dfdt[], void *params);
    double y_exact(double x, double y0);

    double x0 = 1.0;
    double x1 = 10.;
    double step = 1e-3;
    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    std::vector<double> init_values = {0.0, 1.0, 5.0, 7.5};
};

class Solver_ex1 {
    
    private:
    static int func(double x, const double y[], double f[], void *params);
    
    double t0 = 0.0;
    int k = 1;
    double t1 = 3.;
    double step = 1e-3;
    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    // std::vector<double> init_values = {-1., -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    std::vector<double> init_values = {-0.1, 0.1, 0.9, 1.1};
    public:
        Solver_ex1();
        void solve_ode();
};

class Solver_ex2 {
    public:
        Solver_ex2();
        void solve_ode(double t0, double t1, int n_steps, const std::string &out_prefix);
        void write_vector_field_2d(double xmin, double xmax, int nx,
                                   double ymin, double ymax, int ny,
                                   const std::string &fname) ;

    private:
        static int func(double x, const double y[], double f[], void *params);
};
