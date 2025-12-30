#pragma once
#include <complex>
#include <vector>
#include <string>
#include <fftw3.h>

struct QuantumDuffingParams {
    double alpha, beta, gamma, omega;  // potential params - without zeta
    double hbar;                       // atomic units
    double mass;                       // atomic units
    double x_min, x_max;               // xlims
    int N;                             // number of grid points
};

class QuantumDuffing {
private:
    QuantumDuffingParams p;
    
    std::vector<double> x; // position
    std::vector<double> k; // momentum space
    double dx;             // step in position space
    double dk;             // step in momentum space
    
    std::vector<std::complex<double>> psi;      // wavefunction
    
    /**
     * @brief FFTW plans
     */
    fftw_plan plan_forward;
    fftw_plan plan_backward;
    fftw_complex *fft_in, *fft_out;

    std::vector<double> V; // potential 
    
    void setup_grid(); // resize the x and k grids
    void setup_fftw(); // alloc memory for fft plans
    void cleanup_fftw(); // destructor
    double potential(double x_val, double t) const; // duffing potential for a given x and t
    void apply_potential_half_step(double dt, double t); // strang half step
    void apply_kinetic_step(double dt);// strang full step
    
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
    
    void set_initial_gaussian(double x0, double p0, double sigma); // initial condition is a gaussian
    /**
     * @brief main function of the evolution
     * - saves data, applies potential half step, kinetic step and potential half step
     */
    void evolve(double dt, int n_steps, 
                std::string filename = "quantum_evolution",
                int save_every = 10);
    // expectation values
    double expectation_x() const;
    double expectation_p() const;
    double expectation_energy(double t) const;
    double norm() const;
    
    void save_state(std::string filename, double t);
};