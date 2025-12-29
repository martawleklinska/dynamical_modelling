#pragma once
#include "moyal_setup.hpp"
#include<memory>

/**
 * @brief gaussian potential
 */
class GaussianPotential : public Potential {
    private:
        double amplitude_;
        double center_;
        double width_;

    public:
        GaussianPotential(double amplitude, double center, double width) : amplitude_(amplitude), center_(center), width_(width) {};

        double operator()(double x, double t = 0.0) const override {
            double dx = x - center_;
            return amplitude_ * std::exp(-dx*dx / (2.0 * width_*width_));
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<GaussianPotential>(amplitude_, center_, width_);
        }
};

/**
 * @brief harmonic oscillator potential
 */
class HarmonicPotential : public Potential {
    private:
        double mass_;
        double omega_;
        double center_;

    public:
        HarmonicPotential(double mass, double omega, double center = 0.0) 
            : mass_(mass), omega_(omega), center_(center) {};
        double operator()(double x, double t = 0.0) const override {
            double dx = x - center_;
            return 0.5 * mass_ * omega_ * omega_ * dx * dx;
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<HarmonicPotential>(mass_, omega_, center_);
        }
};

/**
 * @brief DK potential exp
 */
class OGPotential : public Potential {
    public:
        double operator()(double x, double t = 0.0) const override {
            return 0.008 * std::exp(-std::pow(x + 200.0, 2) / 250.);
        }
        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<OGPotential>();
        }
};
/**
 * @brief Duffing potential
 */
class Duffing : public Potential {
public:
    double alpha = -1e-2;
    double beta  =  1e-4;
    double gamma =  0.5;
    double omega =  1.0e-3;


    double operator()(double x, double t = 0.0) const override {
        return 0.5*alpha*x*x + 0.25*beta*x*x*x*x - x*gamma * std::cos(omega*t);  
    }
    
    std::unique_ptr<Potential> clone() const override {
        return std::make_unique<Duffing>();
    }
};

