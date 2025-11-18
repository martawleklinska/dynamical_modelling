#include "duffing.hpp"
#include <iostream>

int main() {
    Duffing du(0.05, 1.0, 1.0, 0.2, 1.0); // zeta, alpha, beta, gamma, omega
    du.poincare_map(0.25,   // gamma
                    200.0,  // discard_transient (s) -- lub 200*periods
                    200,    // n_periods_sample
                    "duffing"); // prefix

    // bifurkacja w gamma
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