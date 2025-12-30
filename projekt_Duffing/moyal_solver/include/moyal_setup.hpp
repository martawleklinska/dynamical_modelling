#pragma once
#include<vector>
#include<complex>
#include<functional>
#include<memory>

using Complex = std::complex<double>;
using ComplexMatrix = std::vector<std::vector<Complex>>;
using RealMatrix = std::vector<std::vector<double>>;

/**
 * @brief Potential functions class
 * destructor must be virtual
 */
class Potential {
    public:
        virtual ~Potential() = default;
        virtual double operator()(double x, double t = 0.0) const = 0;
        virtual std::unique_ptr<Potential> clone() const = 0;
};

/**
 * @brief config params for Moyal eqn solver
 */

 struct MoyalConfig {
    // grid params
    int gridX = 1024;     //!< position grid points
    int gridP = 512;      //!< momentum grid points
    double ampX = 1500.0; //!< position range
    double ampP = 0.5;    //!< momentum range

    double hbar = 1.0;    //!< hbar - may be reduced if semiclassical is studied
    double mass = 1.0;    //!< electron mass

    // time params
    double dt = 10.;       //!< time step
    int timeSteps = 800;    //!< number of time steps


    // init conditions
    double x_init  = -700.0;
    double p_init  = 0.15;

    double sigma_x = 100.0;
    double sigma_p = 0.005;   


    // output
    int outputEvery = 10; // output frequency
    std::string outputDir = "./output/";
};
