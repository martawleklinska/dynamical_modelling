#pragma once
#include <complex>
#include <vector>
#include <string>
#include <fftw3.h>

struct QuantumDuffingParams {
    double alpha, beta, gamma, omega;  
    double hbar;    
    double mass;    
    double x_min, x_max;
    int N;          
};

class QuantumDuffing {
private:
    QuantumDuffingParams p;
    
    std::vector<double> x;
    std::vector<double> k;
    double dx;
    double dk;
    
    std::vector<std::complex<double>> psi;
    
    fftw_plan plan_forward;
    fftw_plan plan_backward;
    fftw_complex *fft_in, *fft_out;

    std::vector<double> V; 
    
    void setup_grid();
    void setup_fftw();
    void cleanup_fftw();
    double potential(double x_val, double t) const;
    void apply_potential_half_step(double dt, double t);
    void apply_kinetic_step(double dt);
    
public:
    QuantumDuffing(
        double alpha = -1.0,
        double beta = 1.0,
        double gamma = 0.3,
        double omega = 1.0,
        double hbar = 1.0,
        double mass = 1.0,
        double x_min = -4.0,
        double x_max = 4.0,
        int N = 512
    );
    
    ~QuantumDuffing();
    
    void set_initial_gaussian(double x0, double p0, double sigma);
    
    void evolve(double dt, int n_steps, 
                std::string filename = "quantum_evolution",
                int save_every = 10);
    
    double expectation_x() const;
    double expectation_p() const;
    double expectation_energy(double t) const;
    double norm() const;
    
    void save_state(std::string filename, double t);
};