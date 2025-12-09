#include"hopf.hpp"
#include<iostream>
#include<vector>
#include<array>

int main(){

    for (int i = 0; i<5; i++){
        double m = 0.5 + i * 0.5;
        SolverVdP ex1(0.5, m);
            std::vector<std::array<double,2>> initial_conditions = {
                {0.1, 0.2},   // x0=0.1, y0=0.2
                {-0.5, 1.0},  // x0=-0.5, y0=1.0
                {2.0, -1.5},  // x0=2.0, y0=-1.5
                {0.0, 0.0}    // x0=0.0, y0=0.0
            };

            ex1.solve_trajectories(initial_conditions, "eps05_m"+std::to_string(m));}

    return 0;
}