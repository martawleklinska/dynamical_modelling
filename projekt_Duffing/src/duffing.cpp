#include "duffing.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <array>

Duffing::Duffing(double zeta_, double alpha_, double beta_, double gamma_, double omega_) {
    p = {zeta_, alpha_, beta_, gamma_, omega_};
}

int Duffing::func(double t, const double y[], double f[], void *params) {
    auto *p = static_cast<DuffingParams*>(params);
    double x = y[0];
    double v = y[1];

    f[0] = v;
    f[1] = -2.0 * p->zeta * v - p->alpha * x - p->beta * std::pow(x, 3) + p->gamma * std::cos(p->omega * t);
    return GSL_SUCCESS;
}

void Duffing::solve(std::string filename) {
    std::filesystem::create_directories("data");

    gsl_odeiv2_system sys = {func, nullptr, 2, &p};

    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, h, 1e-8, 1e-8);


    for (auto ic : initials) {
        std::string fname = "data/duff_" + filename +  "_x0_" + std::to_string(ic[0]) +
                            "_v0_" + std::to_string(ic[1]) + ".txt";
        std::ofstream ofs(fname);
        ofs << "t\tx\tv\n";

        double y[2] = {ic[0], ic[1]};
        double t = t0;

        for (int i = 0; i <= n_steps; ++i) {
            double ti = t0 + i * (t1 - t0) / n_steps;
            int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
            if (status != GSL_SUCCESS) {
                std::cerr << "GSL error at step " << i << "\n";
                break;
            }
            ofs << ti << "\t" << y[0] << "\t" << y[1] << "\n";
        }
        ofs.close();
    }

    gsl_odeiv2_driver_free(d);
}

void Duffing::poincare_map(double discard_transient, int n_periods_sample, const std::string &out_prefix, double t0_override){
    double T = 2.0 * M_PI / p.omega;

    gsl_odeiv2_system sys = {func, nullptr, 2, &p};
    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-09, 1e-09);
    
    std::vector<std::array<double, 2>> init_vals = {
        {0.5, 0.0}
    };
    std::string fname = "data/" + out_prefix + "_poincare.txt";
    std::ofstream ofs(fname);
    ofs << "period_idx\t x\t v \n";
    for (auto ic : init_vals){
        double y[2] = {ic[0], ic[1]};
        double t = (t0_override >= 0.0 ? t0_override : t0);

        double t_end_trans = t+ discard_transient;
        int status = gsl_odeiv2_driver_apply(d, &t, t_end_trans, y);
        if (status != GSL_SUCCESS){
            std::cerr << "GSL error during transient\n";
            gsl_odeiv2_driver_free(d);
            return ;
        }
        for (int n = 1; n <= n_periods_sample; ++n){
            double t_target = t_end_trans + n * T;
            status = gsl_odeiv2_driver_apply(d, &t, t_target, y);
            if (status != GSL_SUCCESS){
                std::cerr << "GLS error furing sampling\n";
                break ;
            }
            ofs << n << "\t" << y[0] << "\t" << y[1] << "\n"; 
        }
    }
    ofs.close();
    gsl_odeiv2_driver_free(d);
}


void Duffing::bifurcation_analysis(
    std::string param_name,
    double param_min,
    double param_max,
    int n_param_steps,
    std::array<double, 2> ic,
    double t_transient,
    int n_periods,
    std::string filename
) {
    std::filesystem::create_directories("data");
    
    std::string fname = "data/bifurcation_" + filename + "_" + param_name + ".txt";
    std::ofstream ofs(fname);
    ofs << "# Bifurcation analysis: parameter " << param_name << "\n";
    ofs << "# param\tx\tv\n";
    
    DuffingParams p_original = p;
    
    double param_step = (param_max - param_min) / (n_param_steps - 1);
    
    for (int ip = 0; ip < n_param_steps; ++ip) {
        double param_value = param_min + ip * param_step;
        
        if (param_name == "zeta") p.zeta = param_value;
        else if (param_name == "alpha") p.alpha = param_value;
        else if (param_name == "beta") p.beta = param_value;
        else if (param_name == "gamma") p.gamma = param_value;
        else if (param_name == "omega") p.omega = param_value;
        else {
            std::cerr << "Unknown parameter: " << param_name << "\n";
            return;
        }
        
        gsl_odeiv2_system sys = {func, nullptr, 2, &p};
        double h = 2.0 * M_PI / p.omega / 100.0;  
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
            &sys, gsl_odeiv2_step_rk8pd, h, 1e-10, 1e-10);
        
        double y[2] = {ic[0], ic[1]};
        double t = 0.0;
        
        int status = gsl_odeiv2_driver_apply(d, &t, t_transient, y);
        if (status != GSL_SUCCESS) {
            std::cerr << "GSL error during transient at param=" << param_value << "\n";
            gsl_odeiv2_driver_free(d);
            continue;
        }
        
        double T = 2.0 * M_PI / p.omega;
        
        for (int period = 0; period < n_periods; ++period) {
            double t_next = t + T;
            status = gsl_odeiv2_driver_apply(d, &t, t_next, y);
            
            if (status != GSL_SUCCESS) {
                std::cerr << "GSL error at period " << period 
                          << ", param=" << param_value << "\n";
                break;
            }
            
            ofs << param_value << "\t" << y[0] << "\t" << y[1] << "\n";
        }
        
        gsl_odeiv2_driver_free(d);
        
        if ((ip + 1) % (n_param_steps / 10) == 0 || ip == 0) {
            std::cout << "Progress: " << (ip + 1) << "/" << n_param_steps 
                      << " (" << param_name << " = " << param_value << ")\n";
        }
    }
    
    ofs.close();
    
    p = p_original;
    
    std::cout << "Bifurcation analysis saved to: " << fname << "\n";
}

double Duffing::lyapunov_exponent(
    std::array<double, 2> ic,
    double d0,
    double t_transient,
    double t_measure,
    int n_renorm
) {
    gsl_odeiv2_system sys = {func, nullptr, 2, &p};
    double h = 2.0 * M_PI / p.omega / 100.0;
    
    gsl_odeiv2_driver *d1 = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, h, 1e-10, 1e-10);
    gsl_odeiv2_driver *d2 = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, h, 1e-10, 1e-10);
    
    double y1[2] = {ic[0], ic[1]};
    double y2[2] = {ic[0] + d0, ic[1]};
    
    double t1 = 0.0, t2 = 0.0;
    
    gsl_odeiv2_driver_apply(d1, &t1, t_transient, y1);
    gsl_odeiv2_driver_apply(d2, &t2, t_transient, y2);
    
    double sum_log = 0.0;
    double dt = t_measure / n_renorm;
    
    for (int i = 0; i < n_renorm; ++i) {
        double t_next = t1 + dt;
        
        gsl_odeiv2_driver_apply(d1, &t1, t_next, y1);
        gsl_odeiv2_driver_apply(d2, &t2, t_next, y2);
        
        double dx = y2[0] - y1[0];
        double dv = y2[1] - y1[1];
        double d = std::sqrt(dx*dx + dv*dv);
        
        sum_log += std::log(d / d0);
        
        double scale = d0 / d;
        y2[0] = y1[0] + dx * scale;
        y2[1] = y1[1] + dv * scale;
    }
    
    gsl_odeiv2_driver_free(d1);
    gsl_odeiv2_driver_free(d2);
    
    return sum_log / t_measure;
}

void Duffing::lyapunov_vs_parameter(
    std::string param_name,
    double param_min,
    double param_max,
    int n_steps,
    std::array<double, 2> ic,
    std::string filename
) {
    std::filesystem::create_directories("data");
    
    std::string fname = "data/lyapunov_" + filename + "_" + param_name + ".txt";
    std::ofstream ofs(fname);
    ofs << "# Lyapunov exponent vs " << param_name << "\n";
    ofs << "# param\tlambda\n";
    
    DuffingParams p_original = p;
    double param_step = (param_max - param_min) / (n_steps - 1);
    
    for (int i = 0; i < n_steps; ++i) {
        double param_value = param_min + i * param_step;
        
        if (param_name == "zeta") p.zeta = param_value;
        else if (param_name == "alpha") p.alpha = param_value;
        else if (param_name == "beta") p.beta = param_value;
        else if (param_name == "gamma") p.gamma = param_value;
        else if (param_name == "omega") p.omega = param_value;
        
        double lambda = lyapunov_exponent(ic, 1e-8, 200.0, 500.0, 50000);
        
        ofs << param_value << "\t" << lambda << "\n";
        
        if ((i + 1) % (n_steps / 10) == 0 || i == 0) {
            std::cout << "Progress: " << (i + 1) << "/" << n_steps 
                      << " (" << param_name << " = " << param_value 
                      << ", λ = " << lambda << ")\n";
        }
    }
    
    ofs.close();
    p = p_original;
    
    std::cout << "Lyapunov analysis saved to: " << fname << "\n";
}


static int func_with_energy(double t, const double y[], double f[], void *params) {
    auto *p = static_cast<DuffingParamsWithEnergy*>(params);
    double x = y[0];
    double v = y[1];
    double E = y[2];

    f[0] = v;
    f[1] = -2.0 * p->zeta * v - p->alpha * x - p->beta * std::pow(x, 3) 
           + p->gamma * std::cos(p->omega * t);
    
    f[2] = p->gamma * std::cos(p->omega * t) * v;
    
    return GSL_SUCCESS;
}

void Duffing::solve_with_energy(std::string filename, std::array<double, 2> ic) {
    std::filesystem::create_directories("data");
    
    std::string fname = "data/energy_" + filename + ".txt";
    std::ofstream ofs(fname);
    ofs << "# t\tx\tv\tE_kinetic\tE_potential\tE_input\tE_total\n";
    
    DuffingParamsWithEnergy p_with_E = {p.zeta, p.alpha, p.beta, p.gamma, p.omega, nullptr};
    
    gsl_odeiv2_system sys = {func_with_energy, nullptr, 3, &p_with_E};
    
    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, h, 1e-10, 1e-10);
    
    double y[3] = {ic[0], ic[1], 0.0};
    double t = t0;
    
    for (int i = 0; i <= n_steps; ++i) {
        double ti = t0 + i * (t1 - t0) / n_steps;
        gsl_odeiv2_driver_apply(d, &t, ti, y);
        
        double x = y[0];
        double v = y[1];
        double E_input = y[2];
        
        double E_kinetic = 0.5 * v * v;
        double E_potential = 0.5 * p.alpha * x * x + 0.25 * p.beta * x * x * x * x;
        double E_total = E_kinetic + E_potential;
        
        ofs << ti << "\t" << x << "\t" << v << "\t" 
            << E_kinetic << "\t" << E_potential << "\t" 
            << E_input << "\t" << E_total << "\n";
    }
    
    ofs.close();
    gsl_odeiv2_driver_free(d);
    
    std::cout << "Energy analysis saved to: " << fname << "\n";
}