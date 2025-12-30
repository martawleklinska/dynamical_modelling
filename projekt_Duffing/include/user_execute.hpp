#ifndef USER_EXECUTE_HPP
#define USER_EXECUTE_HPP

/**
 * @brief functions that initializes the calculations for a given analysis
 */

/** 
 * @brief with gamma=0, omega=0 the function solves the duffing oscillator
 * for few most important cases of zeta, alpha, beta configurations - data saved to build/data
 */
void get_time_dependent_solver();

/** 
 * @brief get data to plot Poincare maps (already the params are set to give chaos)
 */
void run_poincare_analysis();
/** 
 * @brief get data to plot bifurcation diagram or gif the Poincare maps
 */
void run_bifurcation_analysis();
/** 
 * @brief get data to plot Lapunow exponent in gamma (to get other params change it in .cpp)
 */
void run_lyapunov_analysis();
/** 
 * @brief get data to plot energy analysis (setails of the params in .cpp)
 */
void run_energy_analysis();
/** 
 * @brief time dependent Schrodinger equation solver
 */
void run_quantum_analysis();
/** 
 * @brief show the user gui menu in the terminal
 */
void show_menu_and_execute();

#endif
