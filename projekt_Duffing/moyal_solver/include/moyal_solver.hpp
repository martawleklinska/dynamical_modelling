#pragma once
#include"moyal_setup.hpp"
#include"phase_space.hpp"
#include"propagators.hpp"
#include"wigner_df.hpp"

/**
 * @brief main Moyal equation solver using Strang splittling
 */

 class MoyalSolver{
    private:
        MoyalConfig config_;
        PhaseSpace phase_space_;
        WDF wigner_;

        std::unique_ptr<KineticPropagator> kinetic_prop_;
        std::unique_ptr<PotentialPropagator> potential_prop_;

        int current_step_;
        double current_time_;

    public:
        explicit MoyalSolver(const MoyalConfig& config, std::unique_ptr<Potential> potential);

        // main evolution methods
        void setInitialCondition(const std::function<Complex(double, double)>& init_func);
        void setInitialGaussian(double x0, double p0, double sigma_x, double sigma_p);

        void evolve(int steps = -1); // -1 means config_.timeSteps
        void evolveOneStep();

        //access
        const WDF& getWigner() const { return wigner_ ;}
        const PhaseSpace& getPhaseSpace() const { return phase_space_;}
        double getCurrentTime() const { return current_time_; }
        int getCurrentStep() const { return current_step_;}

        //analysis
        void computeExpectationValues(double& mean_x, double& mean_p,
                                        double& sigma_x, double& sigma_p) const;
        double computeEnergy() const;

        // saving
        void saveState(const std::string& filename) const;
        void loadState(const std::string& filename);
        void setupOutput();

    private:
        void strangSplittingStep();
        void outputData();
 }; 