#pragma once
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <string>
#include<vector>
#include<array>

/**
 * @brief Duffing Params struct 
 */
struct DuffingParams {
    double zeta, alpha, beta, gamma, omega;
};

class Duffing {
private:
    // function in the Right Hand Side of the differential equation of 1st order
    static int func(double t, const double y[], double f[], void *params);
    DuffingParams p; 

    double t0 = 0.0; // start 
    double t1 = 10.; // end
    int n_steps = 10000; // n/o of steps
    std::vector<std::array<double, 2>> initials = {         // initial conditions of (x0, v0)=(x0, \dot{x}0)
        {0.5, -2.0}, {-1.0, 2.0}, {-1., 0.50}, {0.5, -1.70},
        {0.0, 0.20}, {0.0, 0.50}};

public:
    /**
     * @brief Constructor - default Duffing params
     */
    Duffing(double zeta = .05, double alpha = 1., double beta = 1., 
            double gamma = .2, double omega = 1.);

    /**
     * @brief main solver - using GSL library solving the DE with set Duffing Params
     */
    void solve(std::string filename);

    /**
     * @brief solve and save data to poincare map (gamma)
     */
    void poincare_map(double discard_transient,
                      int n_periods_sample,
                      const std::string &out_prefix,
                      double t0_override = -1.0);

    /**
     * @brief bifurcation analysis - it overlaps with poincare 
     * we can change different params
     * @param param_name - we can chose what params to change
     * @param param_min, param_max - limits
     */
    void bifurcation_analysis(
        std::string param_name,      // name of the parameter: "zeta", "alpha", "beta", "gamma", "omega"
        double param_min,             // minimal limit of param
        double param_max,             // max limit of the param
        int n_param_steps,            // number of param steps
        std::array<double, 2> ic,     // initial condition (x0, v0)
        double t_transient = 100.0,   // transient time
        int n_periods = 50,           // number of periods to save (poincare)
        std::string filename = "bifurcation"
    );
    /**
     * @brief Lapunow exponent analysis
     */
    double lyapunov_exponent(
        std::array<double, 2> ic,     // initial condition
        double d0 = 1e-8,              // difference of initial condition 
        double t_transient = 100.0,    // transient time
        double t_measure = 500.0,      // time of measure
        int n_renorm = 50000           // renormalisation factor
    );
    /**
     * @brief using lyapunov_exponent() we can check its dependency on a given param
     */
    void lyapunov_vs_parameter(
        std::string param_name,// name of the parameter: "zeta", "alpha", "beta", "gamma", "omega"
        double param_min,
        double param_max,
        int n_steps,
        std::array<double, 2> ic,
        std::string filename = "lyapunov"
    );
    /**
     * @brief energy analysis
     */
    void solve_with_energy(
        std::string filename,
        std::array<double, 2> ic = {0.0, 0.0}
    );

    void resonance_curve(
        double omega_min,
        double omega_max,
        int n_steps,
        std::array<double, 2> ic,
        double t_transient = 200.0,
        int n_periods_measure = 50,
        std::string filename = "resonance"
    );
};

struct DuffingParamsWithEnergy {
    double zeta, alpha, beta, gamma, omega;
    double *energy_ptr;  
};