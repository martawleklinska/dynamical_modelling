#pragma once
#include <vector>
#include <string>
#include <complex>
#include <fftw3.h>

struct WignerParams {
    double alpha, beta, gamma, omega;  
    double hbar, mass;                 
    double x_min, x_max, p_min, p_max;
    int Nx, Np;
};

class WignerDuffing {
private:
    WignerParams p;
    
    std::vector<double> x, p_vals;  
    std::vector<double> lambda, theta;
    double dx, dp;
    
    std::vector<std::vector<std::complex<double>>> W;
    
    fftw_complex *data_in, *data_out;
    fftw_plan plan_x_fwd, plan_x_bwd;
    fftw_plan plan_p_fwd, plan_p_bwd;
    
    void setup_grid();
    void setup_fftw();
    void cleanup_fftw();
    
    double U_operator(double x_val, double theta_val, double t) const;
    
    void fft_in_x_direction(bool forward);
    void fft_in_p_direction(bool forward);
    
    void step_kinetic(double dt);
    void step_potential(double dt, double t);
    
public:
    WignerDuffing(double alpha, double beta, double gamma, double omega,
                  double hbar, double mass,
                  double x_min, double x_max, 
                  double p_min, double p_max,
                  int Nx, int Np);
    ~WignerDuffing();
    
    void set_initial_gaussian(double x0, double p0, double sigma);
    void evolve(double dt, int n_steps, std::string filename, int save_every);
    
    // Obserwable
    double expectation_x() const;
    double expectation_p() const;
    double total_probability() const;
    
    void save_state(std::string filename, double t);
};