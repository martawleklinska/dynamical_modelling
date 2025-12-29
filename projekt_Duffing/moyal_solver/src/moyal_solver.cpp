#include "moyal_solver.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

MoyalSolver::MoyalSolver(const MoyalConfig& config, std::unique_ptr<Potential> potential)
    : config_(config), phase_space_(config), wigner_(phase_space_),
      current_step_(0), current_time_(0.0) {
    
    kinetic_prop_ = std::make_unique<KineticPropagator>(phase_space_, config_.mass);
    potential_prop_ = std::make_unique<PotentialPropagator>(phase_space_, 
                                                           std::move(potential), 
                                                           config_.hbar);
}

void MoyalSolver::setInitialCondition(const std::function<Complex(double, double)>& init_func) {
    wigner_.initializeFromFunction(init_func);
}

void MoyalSolver::setInitialGaussian(double x0, double p0, double sigma_x, double sigma_p) {
    wigner_.initializeGaussian(x0, p0, sigma_x, sigma_p);
}

void MoyalSolver::evolve(int steps) {
    int total_steps = (steps > 0) ? steps : config_.timeSteps;
    
    std::cout << "Starting evolution for " << total_steps << " steps...\n";
    
    // Save initial state
    if (current_step_ % config_.outputEvery == 0) {
        outputData();
    }
    
    for (int step = 0; step < total_steps; ++step) {
        evolveOneStep();
        
        if ((current_step_) % config_.outputEvery == 0) {
            outputData();
        }
        
        if (current_step_ % 100 == 0) {
            std::cout << "Step " << current_step_ << "/" << total_steps 
                     << ", time = " << current_time_ << std::endl;
        }
    }
    
    std::cout << "Evolution completed.\n";
}

void MoyalSolver::evolveOneStep() {
    strangSplittingStep();
    current_step_++;
    current_time_ += config_.dt;
}

void MoyalSolver::strangSplittingStep() {
    // Strang splitting: exp(-i*dt/2*T) exp(-i*dt*V) exp(-i*dt/2*T)
    
    auto& wigner_data = wigner_.data();
    
    // 1. Half kinetic step
    phase_space_.fft_p(wigner_data, true);  // FFT in momentum direction
    kinetic_prop_->apply(wigner_data, config_.dt);
    phase_space_.fft_p(wigner_data, false); // IFFT
    
    // 2. Full potential step
    phase_space_.fft_x(wigner_data, true);  // FFT in position direction
    potential_prop_->apply(wigner_data, config_.dt);
    phase_space_.fft_x(wigner_data, false); // IFFT
    
    // 3. Half kinetic step
    phase_space_.fft_p(wigner_data, true);  // FFT in momentum direction
    kinetic_prop_->apply(wigner_data, config_.dt);
    phase_space_.fft_p(wigner_data, false); // IFFT
    
    // Keep only real part for Wigner function
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            wigner_data[i][j] = Complex(wigner_data[i][j].real(), 0.0);
        }
    }
}

void MoyalSolver::computeExpectationValues(double& mean_x, double& mean_p, 
                                         double& sigma_x, double& sigma_p) const {
    auto pos_result = wigner_.expectedPosition();
    auto mom_result = wigner_.expectedMomentum();
    
    mean_x = pos_result.first;
    mean_p = mom_result.first;
    
    // Calculate uncertainties (simplified)
    const auto& X = phase_space_.X();
    const auto& P = phase_space_.P();
    const auto& wigner_data = wigner_.data();
    
    int nx = phase_space_.gridX();
    int np = phase_space_.gridP();
    double dx = phase_space_.dx();
    double dp = phase_space_.dp();
    
    double var_x = 0.0, var_p = 0.0, total_prob = 0.0;
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < np; ++j) {
            double w_real = wigner_data[i][j].real();
            if (w_real > 0) {
                double dx_dev = X[i][j] - mean_x;
                double dp_dev = P[i][j] - mean_p;
                
                var_x += dx_dev * dx_dev * w_real * dx * dp;
                var_p += dp_dev * dp_dev * w_real * dx * dp;
                total_prob += w_real * dx * dp;
            }
        }
    }
    
    if (total_prob > 0) {
        sigma_x = std::sqrt(var_x / total_prob);
        sigma_p = std::sqrt(var_p / total_prob);
    } else {
        sigma_x = sigma_p = 0.0;
    }
}

void MoyalSolver::saveState(const std::string& filename) const {
    wigner_.saveToFile(filename);
}

void MoyalSolver::outputData() {
    // Save current Wigner distribution
    std::ostringstream filename;
    filename << config_.outputDir << "wigner_" << std::setw(6) << std::setfill('0') 
             << current_step_ << ".dat";
    
    saveState(filename.str());
    
    // Calculate and save expectation values
    double mean_x, mean_p, sigma_x, sigma_p;
    computeExpectationValues(mean_x, mean_p, sigma_x, sigma_p);
    
    std::ofstream stats_file(config_.outputDir + "stats.dat", std::ios::app);
    if (current_step_ == 0) {
        stats_file << "# Step Time Mean_X Mean_P Sigma_X Sigma_P Norm\n";
    }
    
    stats_file << current_step_ << " " << current_time_ << " "
               << mean_x << " " << mean_p << " " 
               << sigma_x << " " << sigma_p << " "
               << wigner_.totalProbability() << "\n";
}
