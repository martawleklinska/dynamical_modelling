#include"hopf.hpp"
#include<iostream>
#include<vector>
#include<array>

int main(){
    SolverVdP ex1(0.5, 1.0);
    std::vector<std::array<double,2>> initial_conditions = {
        {0.1, 0.2},   // x0=0.1, y0=0.2
        {-0.5, 1.0},  // x0=-0.5, y0=1.0
        {2.0, -1.5},  // x0=2.0, y0=-1.5
        {0.0, 0.0}    // x0=0.0, y0=0.0
    };

    ex1.solve_trajectories(initial_conditions, "eps05m1");
    return 0;
}