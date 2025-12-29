#pragma once
#include "moyal_setup.hpp"
#include "phase_space.hpp"

/**
 * @brief class for evolution propagators
 */
class Propagator {
    protected:
        const PhaseSpace& phase_space_;
        ComplexMatrix propagator_matrix_;
    
    public:
        explicit Propagator(const PhaseSpace& phase_space) : phase_space_(phase_space) {
            int nx = phase_space_.gridX();
            int np = phase_space_.gridP();
            propagator_matrix_.resize(nx, std::vector<Complex>(np));
        }

        virtual ~Propagator() = default;

        // apply propagator to WDF
        virtual void apply(ComplexMatrix& wigner, double dt) = 0;

        // access to propagation matrix
        const ComplexMatrix& matrix() const {return propagator_matrix_;}
    
    protected:
        virtual void computePropagator(double dt) = 0;

};

/**
 * @brief kinetic energy propagator
 */
class KineticPropagator : public Propagator {
    private:
        double mass_;

    public:
        KineticPropagator(const PhaseSpace& phase_space, double mass);

        void apply(ComplexMatrix& wigner, double dt) override;

    protected:
        void computePropagator(double dt) override;
};

/**
 * @brief potential energy propagator with Moyal bracket
 */
class PotentialPropagator : public Propagator {
    private:
        std::unique_ptr<Potential> potential_;
        double hbar_;

    public:
        PotentialPropagator(const PhaseSpace& phase_space,
                            std::unique_ptr<Potential> potential,
                            double hbar = 1.0);
        void apply(ComplexMatrix& wigner, double dt) override;

    protected:
        void computePropagator(double dt) override;

    private:
        void computeMoyalBracket(ComplexMatrix& wigenr, double dt, double t = 0.0);
};