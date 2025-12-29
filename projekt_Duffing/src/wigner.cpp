#include "wigner.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <filesystem>

WignerDuffing::WignerDuffing(
    double alpha, double beta, double gamma, double omega,
    double hbar, double mass,
    double x_min, double x_max, 
    double p_min, double p_max,
    int Nx, int Np
) {
    p = {alpha, beta, gamma, omega, hbar, mass, 
         x_min, x_max, p_min, p_max, Nx, Np};
    
    setup_grid();
    setup_fftw();
    
    // Alokuj W jako kompleksową tablicę
    W.resize(p.Nx, std::vector<std::complex<double>>(p.Np, 0.0));
}

WignerDuffing::~WignerDuffing() {
    cleanup_fftw();
}

void WignerDuffing::setup_grid() {
    dx = (p.x_max - p.x_min) / p.Nx;
    dp = (p.p_max - p.p_min) / p.Np;
    
    x.resize(p.Nx);
    p_vals.resize(p.Np);
    lambda.resize(p.Nx);
    theta.resize(p.Np);
    
    // Siatka x i p (centrowana)
    for (int i = 0; i < p.Nx; ++i) {
        x[i] = p.x_min + (i + 0.5) * dx;
    }
    for (int j = 0; j < p.Np; ++j) {
        p_vals[j] = p.p_min + (j + 0.5) * dp;
    }
    
    // Częstości FFT (fftshift convention)
    double dlambda = 2.0 * M_PI / (p.Nx * dx);
    double dtheta = 2.0 * M_PI / (p.Np * dp);
    
    for (int i = 0; i < p.Nx; ++i) {
        int k = (i < p.Nx/2) ? i : i - p.Nx;
        lambda[i] = k * dlambda;
    }
    for (int j = 0; j < p.Np; ++j) {
        int k = (j < p.Np/2) ? j : j - p.Np;
        theta[j] = k * dtheta;
    }
}

void WignerDuffing::setup_fftw() {
    size_t total_size = p.Nx * p.Np;
    data_in = fftw_alloc_complex(total_size);
    data_out = fftw_alloc_complex(total_size);
    
    // 1D FFT w kierunku x (wiele rzędów)
    int n_x[] = {p.Nx};
    plan_x_fwd = fftw_plan_many_dft(
        1, n_x, p.Np,
        data_in, NULL, 1, p.Nx,
        data_out, NULL, 1, p.Nx,
        FFTW_FORWARD, FFTW_MEASURE
    );
    plan_x_bwd = fftw_plan_many_dft(
        1, n_x, p.Np,
        data_in, NULL, 1, p.Nx,
        data_out, NULL, 1, p.Nx,
        FFTW_BACKWARD, FFTW_MEASURE
    );
    
    // 1D FFT w kierunku p (wiele kolumn)
    int n_p[] = {p.Np};
    plan_p_fwd = fftw_plan_many_dft(
        1, n_p, p.Nx,
        data_in, NULL, p.Nx, 1,
        data_out, NULL, p.Nx, 1,
        FFTW_FORWARD, FFTW_MEASURE
    );
    plan_p_bwd = fftw_plan_many_dft(
        1, n_p, p.Nx,
        data_in, NULL, p.Nx, 1,
        data_out, NULL, p.Nx, 1,
        FFTW_BACKWARD, FFTW_MEASURE
    );
}

void WignerDuffing::cleanup_fftw() {
    fftw_destroy_plan(plan_x_fwd);
    fftw_destroy_plan(plan_x_bwd);
    fftw_destroy_plan(plan_p_fwd);
    fftw_destroy_plan(plan_p_bwd);
    fftw_free(data_in);
    fftw_free(data_out);
}

void WignerDuffing::fft_in_x_direction(bool forward) {
    // Kopiuj W → data_in (row-major: j*Nx + i)
    for (int j = 0; j < p.Np; ++j) {
        for (int i = 0; i < p.Nx; ++i) {
            int idx = j * p.Nx + i;
            data_in[idx][0] = W[i][j].real();
            data_in[idx][1] = W[i][j].imag();
        }
    }
    
    // Wykonaj FFT
    if (forward) {
        fftw_execute(plan_x_fwd);
    } else {
        fftw_execute(plan_x_bwd);
    }
    
    // Kopiuj data_out → W
    double norm = forward ? 1.0 : 1.0 / p.Nx;
    for (int j = 0; j < p.Np; ++j) {
        for (int i = 0; i < p.Nx; ++i) {
            int idx = j * p.Nx + i;
            W[i][j] = norm * std::complex<double>(
                data_out[idx][0], 
                data_out[idx][1]
            );
        }
    }
}

void WignerDuffing::fft_in_p_direction(bool forward) {
    // Kopiuj W → data_in (column-major: i*Np + j)
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            int idx = i * p.Np + j;
            data_in[idx][0] = W[i][j].real();
            data_in[idx][1] = W[i][j].imag();
        }
    }
    
    // Wykonaj FFT
    if (forward) {
        fftw_execute(plan_p_fwd);
    } else {
        fftw_execute(plan_p_bwd);
    }
    
    // Kopiuj data_out → W
    double norm = forward ? 1.0 : 1.0 / p.Np;
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            int idx = i * p.Np + j;
            W[i][j] = norm * std::complex<double>(
                data_out[idx][0], 
                data_out[idx][1]
            );
        }
    }
}

double WignerDuffing::U_operator(double x_val, double theta_val, double t) const {
    // U(x,θ,t) = V(x + ℏθ/2, t) - V(x - ℏθ/2, t)
    
    auto V = [this, t](double x) {
        return 0.5 * p.alpha * x * x 
             + 0.25 * p.beta * x * x * x * x
             - x * p.gamma * std::cos(p.omega * t);
    };
    
    double x_plus = x_val + 0.5 * p.hbar * theta_val;
    double x_minus = x_val - 0.5 * p.hbar * theta_val;
    
    return V(x_plus) - V(x_minus);
}

void WignerDuffing::step_kinetic(double dt) {
    // Mnożenie punktowe w przestrzeni (λ, p):
    // W(λ, p) *= exp[-i dt/(2m) · p²/(2m)]
    //          = exp[-i dt · p²/(4m)]
    
    double coeff = -dt / (4.0 * p.mass);
    
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            double p_val = p_vals[j];
            double phase = coeff * p_val * p_val;
            std::complex<double> factor = std::exp(std::complex<double>(0.0, phase));
            W[i][j] *= factor;
        }
    }
}

void WignerDuffing::step_potential(double dt, double t) {
    // Mnożenie punktowe w przestrzeni (x, θ):
    // W(x, θ) *= exp[-i dt/ℏ · U(x, θ, t)]
    
    double coeff = -dt / p.hbar;
    
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            double U_val = U_operator(x[i], theta[j], t);
            double phase = coeff * U_val;
            std::complex<double> factor = std::exp(std::complex<double>(0.0, phase));
            W[i][j] *= factor;
        }
    }
}

void WignerDuffing::set_initial_gaussian(double x0, double p0, double sigma) {
    double norm = 1.0 / (M_PI * p.hbar);
    
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            double dx_val = x[i] - x0;
            double dp_val = p_vals[j] - p0;
            
            double wigner_val = norm * std::exp(
                -2.0 * dx_val * dx_val / (sigma * sigma)
                -2.0 * sigma * sigma * dp_val * dp_val / (p.hbar * p.hbar)
            );
            
            W[i][j] = std::complex<double>(wigner_val, 0.0);
        }
    }
    
    std::cout << "Initial probability: " << total_probability() << "\n";
}

void WignerDuffing::evolve(double dt, int n_steps, 
                            std::string filename, int save_every) {
    std::filesystem::create_directories("data/wigner");
    
    std::string fname_obs = "data/wigner/" + filename + "_observables.txt";
    std::ofstream ofs_obs(fname_obs);
    ofs_obs << "# step\tt\t<x>\t<p>\tP_total\n";
    
    double t = 0.0;
    
    for (int step = 0; step <= n_steps; ++step) {
        if (step % save_every == 0) {
            std::string fname_state = "data/wigner/" + filename + 
                                       "_step_" + std::to_string(step / save_every) + ".txt";
            save_state(fname_state, t);
            
            double exp_x = expectation_x();
            double exp_p = expectation_p();
            double prob = total_probability();
            
            ofs_obs << step << "\t" << t << "\t" << exp_x << "\t" 
                    << exp_p << "\t" << prob << "\n";
            
            std::cout << "Step " << step << "/" << n_steps 
                      << " t=" << t << ": <x>=" << exp_x 
                      << ", P=" << prob << "\n";
        }
        
        // ===== ALGORYTM 2ND ORDER (Eq. 32) =====
        
        // [1-2] Kinetic half-step
        fft_in_x_direction(true);     // F_{x→λ}: W(x,p) → W(λ,p)
        step_kinetic(dt);              // exp[-i dt p²/(4m)]
        fft_in_x_direction(false);    // F^{λ→x}: W(λ,p) → W(x,p)
        
        // [3-5] Full potential step
        fft_in_p_direction(true);     // F_{p→θ}: W(x,p) → W(x,θ)
        step_potential(dt, t);         // exp[-i dt/ℏ U(x,θ)]
        fft_in_p_direction(false);    // F^{θ→p}: W(x,θ) → W(x,p)
        
        // [6-8] Kinetic half-step
        fft_in_x_direction(true);     // F_{x→λ}
        step_kinetic(dt);              // exp[-i dt p²/(4m)]
        fft_in_x_direction(false);    // F^{λ→x}
        
        t += dt;
    }
    
    ofs_obs.close();
    std::cout << "Wigner evolution completed!\n";
}

double WignerDuffing::total_probability() const {
    double sum = 0.0;
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            sum += W[i][j].real();  // Funkcja Wignera jest rzeczywista
        }
    }
    return sum * dx * dp;
}

double WignerDuffing::expectation_x() const {
    double exp_val = 0.0;
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            exp_val += x[i] * W[i][j].real();
        }
    }
    return exp_val * dx * dp;
}

double WignerDuffing::expectation_p() const {
    double exp_val = 0.0;
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            exp_val += p_vals[j] * W[i][j].real();
        }
    }
    return exp_val * dx * dp;
}

void WignerDuffing::save_state(std::string filename, double t) {
    std::ofstream ofs(filename);
    ofs << "# x\tp\tW(x,p,t=" << t << ")\n";
    
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Np; ++j) {
            ofs << x[i] << "\t" << p_vals[j] << "\t" << W[i][j].real() << "\n";
        }
        ofs << "\n";  // pusta linia dla gnuplot
    }
    
    ofs.close();
}