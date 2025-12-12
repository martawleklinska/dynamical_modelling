#include "duffing.hpp"
#include "user_execute.hpp"
#include <iostream>
#include <cmath>

void get_time_dependent_solver(){
    std::cout << "======= solving the Duffing ODE (γ=0, ω=0) =========\n";
    
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

void show_menu_and_execute() {
    int choice;
    bool continue_running = true;
    
    while (continue_running) {
        std::cout << "\n=========== DUFFING OSCILLATOR MENU ===========\n";
        std::cout << "1. ODE solver (multiple parameter sets)\n";
        std::cout << "2. Poincaré map analysis\n";
        std::cout << "3. Bifurcation analysis\n";
        std::cout << "4. Lyapunov exponent analysis\n";
        std::cout << "0. Exit\n";
        std::cout << "===============================================\n";
        std::cout << "Please enter your choice (0-4): ";
        
        std::cin >> choice;
        
        switch (choice) {
            case 1:
                get_time_dependent_solver();
                break;
            case 2:
                run_poincare_analysis();
                break;
            case 3:
                run_bifurcation_analysis();
                break;
            case 4:
                run_lyapunov_analysis();
                break;
            case 0:
                std::cout << "Exiting... Thank you for using Duffing Oscillator Solver!\n";
                continue_running = false;
                break;
            default:
                std::cout << "Invalid choice! Please enter a number between 0 and 4.\n";
                std::cin.clear();
                std::cin.ignore(10000, '\n');
                break;
        }
        
        if (continue_running && choice >= 1 && choice <= 4) {
            std::cout << "\nPress Enter to return to menu...";
            std::cin.ignore();
            std::cin.get();
        }
    }
}

void run_poincare_analysis() {
    std::cout << "\n======= POINCARÉ MAP ANALYSIS =======\n";
    
    Duffing du(0.1, -1.0, 0.25, 2.5, 2.0); // zeta, alpha, beta, gamma, omega
    
    double T = 2.0 * M_PI / 2.0; 
    std::cout << "Period T = " << T << " seconds" << std::endl;
    
    std::cout << "Running Poincaré map analysis...\n";
    du.poincare_map(10*T,      // discard ~10 periods for transient
                    50000,     // 50000 samples
                    "poincare_analysis");
    
    std::cout << "Poincaré map analysis completed.\n";
}

void run_bifurcation_analysis() {
    std::cout << "\n======= BIFURCATION ANALYSIS =======\n";
    
    Duffing du(0.1, -1.0, 0.25, 2.5, 2.0); // zeta, alpha, beta, gamma, omega
    
    std::cout << "Running bifurcation scan for gamma in [0.0, 1.0]...\n";
    du.bifurcation_scan(0.0, 1.0, 200,   // gamma_min, gamma_max, n_gamma
                        200.0, 50,        // discard_transient, samples_per_gamma
                        "bifurcation_analysis");
    
    std::cout << "Bifurcation analysis completed.\n";
}

void run_lyapunov_analysis() {
    std::cout << "\n======= LYAPUNOV EXPONENT ANALYSIS =======\n";
    
    std::cout << "Parameters: zeta=0.1, alpha=-1.0, beta=0.25, gamma=2.5, omega=2.0\n";
    Duffing duff(0.1, -1.0, 0.25, 2.5, 2.0);
    
    std::cout << "Computing Lyapunov exponents...\n";
    duff.solve("lyapunov_analysis");
    
    std::cout << "Lyapunov exponent analysis completed.\n";
}