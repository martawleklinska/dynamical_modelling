#include "duffing.hpp"
#include "user_execute.hpp"
#include "quantum_duffing.hpp"
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
    {
        std::cout << "zeta=0.05, alpha = -1.0, beta = 0.25\n";
        Duffing duff(0.05, -1.0, 0.25, 2.5, 2.0);
        duff.solve("b_pos_aneg_zeta05_gamma25_omega2");
    }
    {
        std::cout << "zeta=0.05, alpha = 1.0, beta = 5.0\n";
        Duffing duff(0.05, 1.0, 5.0, 2.5, 2.0);
        duff.solve("ab_pos_zeta05_gamma25_omega2");
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
        std::cout << "5. Energy analysis\n";
        std::cout << "6. Quantum Duffing analysis\n";
        std::cout << "0. Exit\n";
        std::cout << "===============================================\n";
        std::cout << "Please enter your choice (0-6): ";
        
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
            case 5:
                run_energy_analysis();
                break;
            case 6:
                run_quantum_analysis();
                break;
            case 0:
                std::cout << "Exiting... Thank you for using Duffing Oscillator Solver!\n";
                continue_running = false;
                break;
            default:
                std::cout << "Invalid choice! Please enter a number between 0 and 6.\n";
                std::cin.clear();
                std::cin.ignore(10000, '\n');
                break;
        }
        
        if (continue_running && choice >= 1 && choice <= 6) {
            std::cout << "\nPress Enter to return to menu...";
            std::cin.ignore();
            std::cin.get();
        }
    }
}

void run_poincare_analysis() {
    std::cout << "\n======= POINCARÉ MAP ANALYSIS =======\n";
    
    // Duffing du(0.1, -1.0, 0.25, 2.5, 2.0); // zeta, alpha, beta, gamma, omega
    Duffing du(0.15, -1.0, 1.0, 0.36, 1.2);
    double T = 2.0 * M_PI / 2.0; 
    std::cout << "Period T = " << T << " seconds" << std::endl;
    
    std::cout << "Running Poincaré map analysis...\n";
    du.poincare_map(10*T,
                    500,
                    "poincare_analysis");
    
    std::cout << "Poincaré map analysis completed.\n";
}

void run_bifurcation_analysis() {
    std::cout << "\n======= BIFURCATION ANALYSIS =======\n";
    
    Duffing du(0.15, -1.0, 1.0, 2.5, 1.2); // zeta, alpha, beta, gamma, omega
    
    std::cout << "Running bifurcation scan for gamma in [0.0, 1.0]...\n";
    du.bifurcation_analysis(
        "gamma",           // parametr do skanowania
        0.25,               // min
        0.5,               // max
        100,              // liczba kroków
        {0.0, 0.0},        // warunek początkowy
        200.0,             // czas transjentowy
        1000,               // liczba okresów do zapisania
        "gamma_scan"       // nazwa pliku
    );
    
    std::cout << "Bifurcation analysis completed.\n";
}

void run_lyapunov_analysis() {
    std::cout << "\n======= LYAPUNOV EXPONENT ANALYSIS =======\n";
    
    Duffing duff(0.15, -1.0, 1.0, 2.5, 1.2);
    duff.lyapunov_vs_parameter("gamma", 0.1, 0.5, 100, {0.1, 0.0});

    std::cout << "\n====== LYAPUNOV EXPONENT ANALYSIS DONE. ========\n";
}

void run_energy_analysis(){
    std::cout << "\n======= ENERGY ANALYSIS =======\n";
    Duffing duff(0.1, -1.0, 5.0, 2.0, 2.0);
    duff.solve_with_energy("energy", {0.1, 0.5});
}

void run_quantum_analysis(){
    std::cout << "\n ======= TDSE ANALYSIS =======\n";
    QuantumDuffing qd(
        -1.0,      // alpha 
        1.0,       // beta
        0.3,       // gamma 
        1.0,       // omega
        1.0,       // hbar 
        1.0,       // mass
        -4.0,      // x_min
        4.0,       // x_max
        512        // N 
    );
    
    qd.set_initial_gaussian(
        1.0,      // x0 
        0.0,       // p0 
        0.3        // sigma 
    );
    
    double dt = 0.005;                    // krok czasowy
    int n_steps = 500;                  // liczba kroków (t_max = 50)
    qd.evolve(dt, n_steps, "duffing_gauss2", 1);  // zapisuj co 1 kroków
    
}
