#pragma once
#include"moyal_setup.hpp"
#include"phase_space.hpp"

/**
 * @brief WDF class
 */
class WDF{
    private:
        const PhaseSpace& phase_space_;
        ComplexMatrix wigner_;

    public:
        explicit WDF(const PhaseSpace& phase_space);

        ComplexMatrix& data() {return wigner_;}
        const ComplexMatrix& data() const {return wigner_;}

        //initialization
        void initializeGaussian(double x0, double p0, double sigma_x, double sigma_p);
        void initializeCoherentState(double x0, double p0, double sigma);
        void initializeFromFunction(const std::function<Complex(double, double)>& func);

        // exp values
        double norm() const;
        double totalProbability() const;
        std::pair<double, double> expectedPosition() const;
        std::pair<double, double> expectedMomentum() const;
        std::pair<double, double> positionUncertainty() const;
        std::pair<double, double> momentumUncertainty() const;
        
        //saving
        void saveToFile(const std::string& filename) const;
        void loadFromFile(const std::string& filename);

        // visualization data
        RealMatrix real() const;
        RealMatrix imag() const;
        RealMatrix magnitude() const;

    private:
        void normalize();
};