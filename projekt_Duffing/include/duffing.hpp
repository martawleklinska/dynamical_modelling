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

    void bifurcation_analysis(
        std::string param_name,      // nazwa parametru: "zeta", "alpha", "beta", "gamma", "omega"
        double param_min,             // minimalna wartość parametru
        double param_max,             // maksymalna wartość parametru
        int n_param_steps,            // liczba kroków parametru
        std::array<double, 2> ic,     // warunek początkowy {x0, v0}
        double t_transient = 100.0,   // czas na ustalenie się (transient)
        int n_periods = 50,           // liczba okresów do zapisania (Poincaré)
        std::string filename = "bifurcation"
    );
    double lyapunov_exponent(
        std::array<double, 2> ic,     // warunek początkowy
        double d0 = 1e-8,              // początkowe zaburzenie
        double t_transient = 100.0,    // czas transjentowy
        double t_measure = 500.0,      // czas pomiaru
        int n_renorm = 50000           // liczba renormalizacji
    );
    
    void lyapunov_vs_parameter(
        std::string param_name,
        double param_min,
        double param_max,
        int n_steps,
        std::array<double, 2> ic,
        std::string filename = "lyapunov"
    );
    void solve_with_energy(
        std::string filename,
        std::array<double, 2> ic = {0.0, 0.0}
    );

};

struct DuffingParamsWithEnergy {
    double zeta, alpha, beta, gamma, omega;
    double *energy_ptr;  
};