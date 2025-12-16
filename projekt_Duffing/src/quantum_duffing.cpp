#include "quantum_duffing.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip>

QuantumDuffing::QuantumDuffing(
    double alpha_, double beta_, double gamma_, double omega_, 
    double hbar_, double mass_, 
    double x_min_, double x_max_, int N_
){
    p = {alpha_, beta_, gamma_, omega_, hbar_, mass_, x_min_, x_max_, N_};

    setup_grid();
    setup_fftw();

    psi.resize(p.N);
    V.resize(p.N);
}


QuantumDuffing::~QuantumDuffing() {
    cleanup_fftw();
}

void QuantumDuffing::setup_grid() {
    dx = (p.x_max - p.x_min) / (p.N - 1);
    dk = 2.0 * M_PI / (p.N * dx);
    
    x.resize(p.N);
    k.resize(p.N);
    
    for (int i = 0; i < p.N; ++i) {
        x[i] = p.x_min + i * dx;
        
        if (i < p.N / 2) {
            k[i] = i * dk;
        } else {
            k[i] = (i - p.N) * dk;
        }
    }
}

void QuantumDuffing::setup_fftw() {

    fft_in = fftw_alloc_complex(p.N);//malloc: dane w pamięci są obok siebie
    fft_out = fftw_alloc_complex(p.N);

    plan_forward = fftw_plan_dft_1d(p.N, fft_in, fft_out, 
                                     FFTW_FORWARD, FFTW_ESTIMATE);// estimate: potestuj kilka transformat (forward -, backward +)
    plan_backward = fftw_plan_dft_1d(p.N, fft_in, fft_out, 
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
}

void QuantumDuffing::cleanup_fftw() {
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(fft_in);
    fftw_free(fft_out);
}

double QuantumDuffing::potential(double x_val, double t) const {
    return 0.5 * p.alpha * x_val * x_val 
           + 0.25 * p.beta * x_val * x_val * x_val * x_val
           - x_val * p.gamma * std::cos(p.omega * t);
}

void QuantumDuffing::set_initial_gaussian(double x0, double p0, double sigma) {
    double norm_factor = std::pow(2.0 * M_PI * sigma * sigma, -0.25);
    
    for (int i = 0; i < p.N; ++i) {
        double dx_val = x[i] - x0;
        double gauss = norm_factor * std::exp(-dx_val * dx_val / (4.0 * sigma * sigma));
        double phase = p0 * x[i] / p.hbar;
        
        psi[i] = gauss * std::exp(std::complex<double>(0.0, phase));
    }
    
    double n = norm();
    for (auto &val : psi) val /= std::sqrt(n);
}

void QuantumDuffing::apply_potential_half_step(double dt, double t) {
    for (int i = 0; i < p.N; ++i) {
        double V_val = potential(x[i], t);
        std::complex<double> phase(0.0, -V_val * dt / (2.0 * p.hbar));
        psi[i] *= std::exp(phase);
    }
}

void QuantumDuffing::apply_kinetic_step(double dt) {
    for (int i = 0; i < p.N; ++i) {
        fft_in[i][0] = psi[i].real();
        fft_in[i][1] = psi[i].imag();
    }
    
    fftw_execute(plan_forward);
    
    for (int i = 0; i < p.N; ++i) {
        double k_val = k[i];
        double E_kin = p.hbar * p.hbar * k_val * k_val / (2.0 * p.mass);
        std::complex<double> phase(0.0, -E_kin * dt / p.hbar);
        std::complex<double> exp_phase = std::exp(phase);
        
        std::complex<double> psi_k(fft_out[i][0], fft_out[i][1]);
        psi_k *= exp_phase;
        
        fft_in[i][0] = psi_k.real();
        fft_in[i][1] = psi_k.imag();
    }
    
    fftw_execute(plan_backward);
    
    for (int i = 0; i < p.N; ++i) {
        psi[i] = std::complex<double>(fft_out[i][0], fft_out[i][1]) / double(p.N);
    }
}
void QuantumDuffing::evolve(double dt, int n_steps, 
                             std::string filename, int save_every) {
    std::filesystem::create_directories("data/quantum");
    
    std::string fname_obs = "data/quantum/" + filename + "_observables.txt";
    std::ofstream ofs_obs(fname_obs);
    ofs_obs << "# step\tt\t<x>\t<p>\t<E>\tnorm\n";  
    
    double t = 0.0;
    
    std::cout << "Initial norm: " << norm() << "\n";
    
    for (int step = 0; step <= n_steps; ++step) {
        if (step % save_every == 0) {
            std::string fname_state = "data/quantum/" + filename + 
                                       "_step_" + std::to_string(step / save_every) + ".txt";
            save_state(fname_state, t);
            
            double exp_x = expectation_x();
            double exp_p = expectation_p();
            double exp_E = expectation_energy(t);
            double n = norm();
            
            ofs_obs << step << "\t" << t << "\t" << exp_x << "\t" << exp_p << "\t"  
                    << exp_E << "\t" << n << "\n";
            
            std::cout << "Step " << step << "/" << n_steps 
                      << " (t=" << t << "): <x>=" << exp_x 
                      << ", norm=" << n << "\n";
            
            if (std::abs(n - 1.0) > 0.01) {
                std::cerr << "WARNING: Norm deviation at step " << step 
                          << ": " << n << "\n";
            }
        }
        
        apply_potential_half_step(dt, t);
        apply_kinetic_step(dt);
        apply_potential_half_step(dt, t + dt);
        
        t += dt;
    }
    
    ofs_obs.close();
    std::cout << "Quantum evolution completed. Final norm: " << norm() << "\n";
}

double QuantumDuffing::expectation_x() const {
    double exp_val = 0.0;
    for (int i = 0; i < p.N; ++i) {
        exp_val += x[i] * std::norm(psi[i]) * dx;
    }
    return exp_val;
}

double QuantumDuffing::expectation_p() const {
    double exp_val = 0.0;
    for (int i = 1; i < p.N - 1; ++i) {
        std::complex<double> dpsi_dx = (psi[i+1] - psi[i-1]) / (2.0 * dx);
        std::complex<double> integrand = std::conj(psi[i]) * dpsi_dx;
        exp_val += integrand.imag() * dx;
    }
    return p.hbar * exp_val;
}

double QuantumDuffing::expectation_energy(double t) const {
    double E_kin = 0.0;
    double E_pot = 0.0;
    
    for (int i = 1; i < p.N - 1; ++i) {
        std::complex<double> d2psi_dx2 = (psi[i+1] - 2.0*psi[i] + psi[i-1]) / (dx * dx);
        std::complex<double> integrand = std::conj(psi[i]) * d2psi_dx2;
        E_kin += -integrand.real() * dx;
    }
    E_kin *= p.hbar * p.hbar / (2.0 * p.mass);
    
    for (int i = 0; i < p.N; ++i) {
        double V_val = potential(x[i], t);
        E_pot += V_val * std::norm(psi[i]) * dx;
    }
    
    return E_kin + E_pot;
}

double QuantumDuffing::norm() const {
    double n = 0.0;
    for (int i = 0; i < p.N; ++i) {
        n += std::norm(psi[i]) * dx;
    }
    return n;
}

void QuantumDuffing::save_state(std::string filename, double t) {
    std::ofstream ofs(filename);
    ofs << "# x\t|psi|^2\tRe(psi)\tIm(psi)\tV(x,t=" << t << ")\n";
    
    for (int i = 0; i < p.N; ++i) {
        double prob_density = std::norm(psi[i]);
        double V_val = potential(x[i], t);
        
        ofs << x[i] << "\t" << prob_density << "\t" 
            << psi[i].real() << "\t" << psi[i].imag() << "\t"
            << V_val << "\n";
    }
    
    ofs.close();
}