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

    std::vector<std::array<double, 2>> initials = {
        {0.5, -2.03}, {-1.0, 2.03}, {-1., 0.503}, {0.5, -1.703},
        {0.0, 0.203}, {0.0, 0.503}};

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

void Duffing::poincare_map(double discard_transient,
                           int n_periods_sample,
                           const std::string &out_prefix,
                           double t0_override)
{
    double T = 2.0 * M_PI / p.omega;

    gsl_odeiv2_system sys = {func, nullptr, 2, &p};
    double h = (t1 - t0) / n_steps;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-9, 1e-9);

    std::vector<std::array<double,2>> initials = {
        {0.5, 0.0} 
    };

    std::string fname = "data/" + out_prefix + "_poincare.txt";
    std::ofstream ofs(fname);
    ofs << "# period_index\t x\t v\n";
    ofs << std::setprecision(12);

    for (auto ic : initials) {
        double y[2] = {ic[0], ic[1]};
        double t = (t0_override >= 0.0 ? t0_override : t0);

        double t_end_trans = t + discard_transient;
        int status = gsl_odeiv2_driver_apply(d, &t, t_end_trans, y);
        if (status != GSL_SUCCESS) {
            std::cerr << "GSL error during transient (poincare_map) \n";
            gsl_odeiv2_driver_free(d);
            return;
        }

        for (int n = 1; n <= n_periods_sample; ++n) {
            double t_target = t_end_trans + n * T;
            status = gsl_odeiv2_driver_apply(d, &t, t_target, y);
            if (status != GSL_SUCCESS) {
                std::cerr << "GSL error during sampling (poincare_map) \n";
                break;
            }
            ofs << n << "\t" << y[0] << "\t" << y[1] << "\n";
        }
    }

    ofs.close();
    gsl_odeiv2_driver_free(d);
}


void Duffing::bifurcation_scan(double gamma_min, double gamma_max, int n_gamma,
                               double discard_transient, int samples_per_gamma,
                               const std::string &out_prefix,
                               double t0_override)
{
    static const int SUBSTEPS_PER_PERIOD = 200;
    std::string fname_poincare = "data/" + out_prefix + "_bifur_poincare.txt";
    std::string fname_amp = "data/" + out_prefix + "_bifur_amp.txt";

    std::ofstream ofs_p(fname_poincare);
    std::ofstream ofs_a(fname_amp);
    ofs_p << "# gamma\tperiod_index\tx\tv\n";
    ofs_a << "# gamma\tamplitude_max_over_sampled_periods\n";
    ofs_p << std::setprecision(12);
    ofs_a << std::setprecision(12);

    for (int ig = 0; ig < n_gamma; ++ig) {
        double gamma_val = gamma_min + (gamma_max - gamma_min) * ig / double(std::max(1, n_gamma - 1));
        p.gamma = gamma_val;
        double T = 2.0 * M_PI / p.omega;

        gsl_odeiv2_system sys = {func, nullptr, 2, &p};
        double h = (t1 - t0) / n_steps;
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-9, 1e-9);

        double y[2] = {0.5, 0.0};
        double t = (t0_override >= 0.0 ? t0_override : t0);

        double t_end_trans = t + discard_transient;
        int status = gsl_odeiv2_driver_apply(d, &t, t_end_trans, y);
        if (status != GSL_SUCCESS) {
            std::cerr << "GSL error (transient) gamma=" << gamma_val << "\n";
            gsl_odeiv2_driver_free(d);
            continue;
        }

        double amplitude_sum = 0.0;
        int amplitude_count = 0;
        for (int n = 1; n <= samples_per_gamma; ++n) {
            double period_max = -1e300;
            for (int k = 1; k <= SUBSTEPS_PER_PERIOD; ++k) {
                double frac = double(k) / double(SUBSTEPS_PER_PERIOD);
                double t_target = t_end_trans + (n - 1 + frac) * T;
                status = gsl_odeiv2_driver_apply(d, &t, t_target, y);
                if (status != GSL_SUCCESS) {
                    std::cerr << "GSL error (sampling) gamma=" << gamma_val << "\n";
                    break;
                }
                if (y[0] > period_max) period_max = y[0];
            }
            if (period_max > -1e200) {
                amplitude_sum += period_max;
                amplitude_count += 1;
            }

            double t_strobe = t_end_trans + n * T;
            status = gsl_odeiv2_driver_apply(d, &t, t_strobe, y);
            if (status != GSL_SUCCESS) {
                std::cerr << "GSL error (strobe) gamma=" << gamma_val << "\n";
                break;
            }
            ofs_p << gamma_val << "\t" << n << "\t" << y[0] << "\t" << y[1] << "\n";
        }

        double amplitude_avg = (amplitude_count > 0) ? amplitude_sum / amplitude_count : 0.0;
        ofs_a << gamma_val << "\t" << amplitude_avg << "\n";

        gsl_odeiv2_driver_free(d);
    }

    ofs_p.close();
    ofs_a.close();
}