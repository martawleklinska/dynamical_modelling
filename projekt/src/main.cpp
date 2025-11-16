#include "duffing.hpp"
#include <iostream>

int main() {
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
