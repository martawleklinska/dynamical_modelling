#include "duffing.hpp"
#include <iostream>
#include <cmath>

void get_time_dependent_solver();

int main() {
    //
    get_time_dependent_solver();
    //
    // std::cout << "get poincare"<< std::endl;
    // Duffing du(0.1, -1.0, 0.25, 2.5, 2.0); // zeta, alpha, beta, gamma, omega
    
    // double T = 2.0 * M_PI / 2.0;  
    // std::cout << "Period T = " << T << " seconds" << std::endl;

    // du.poincare_map(10*T,   // discard ~10 periods (≈126s) instead of 300s
    //                 50000,     // 50 samples instead of 200
    //                 "poincare1");
    // {
        // std::cout << "zeta=0.1, alpha = 1.0, beta = 5.0\n";
        // Duffing duff(0.1, -1.0, 0.25, 2.5, 2.0);
        // duff.solve("lapunow2");
    // }
    // du.bifurcation_scan(0.0, 1.0, 200,   // gamma_min, gamma_max, n_gamma
    //                     200.0, 50,       // discard_transient, samples_per_gamma
    //                     "duffing");
    return 0;
}

void get_time_dependent_solver(){
    std::cout << "======= solving the Duffing ODE =========\n";
    
    {
        std::cout << "zeta=0.1, alpha = 1.0, beta = 5.0\n";
        Duffing duff(.1, 1.0, 5.0, 0.0, 0.0);
        duff.solve("ab_pos_zeta_small");
    }
    
    {
        std::cout << "zeta=0.1, alpha = -1.0, beta = 5.0\n";
        Duffing duff(.1, -1.0, 5.0, 0.0, 0.0);
        duff.solve("b_pos_a_neg_zeta_small");
    }
    
    {
        std::cout << "zeta=0.1, alpha = 1.0, beta = -5.0\n";
        Duffing duff(.1, 1.0, -5.0, 0.0, 0.0);
        duff.solve("a_pos_b_neg_zeta_small");
    }
    
    {
        std::cout << "zeta=0.1, alpha = -1.0, beta = -5.0\n";
        Duffing duff(.1, -1.0, -5.0, 0.0, 0.0);
        duff.solve("ab_neg_zeta_small");
    }
    
    {
        std::cout << "zeta=3., alpha = -1.0, beta = -5.0\n";
        Duffing duff(3., -1.0, -5.0, 0.0, 0.0);
        duff.solve("ab_neg_zeta_big");
    }
    
    {
        std::cout << "zeta=3., alpha = 1.0, beta = 5.0\n";
        Duffing duff(3., 1.0, 5.0, 0.0, 0.0);
        duff.solve("ab_pos_zeta_big");
    }
    
    {
        std::cout << "zeta=-.1, alpha = 1.0, beta = 5.0\n";
        Duffing duff(-0.1, 1.0, 5.0, 0.0, 0.0);
        duff.solve("ab_pos_zeta_neg");
    }
    
    {
        std::cout << "zeta=-.1, alpha = -1.0, beta = -5.0\n";
        Duffing duff(-0.1, -1.0, -5.0, 0.0, 0.0);
        duff.solve("ab_neg_zeta_neg");
    }
    
    {
        std::cout << "zeta=-3, alpha = -1.0, beta = -5.0\n";
        Duffing duff(-3., -1.0, -5.0, 0.0, 0.0);
        duff.solve("ab_neg_zeta_neg_big");
    }
    
    std::cout << "finished.\n";
}