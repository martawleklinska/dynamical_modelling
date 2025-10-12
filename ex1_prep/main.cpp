// Implement the Euler and Runge–Kutta 4 (RK4) methods for

#include<iostream>
#include<fstream>
#include<cmath>
#include "differentiate.hpp"

int main(){
    double t0 = 0.0;
    double y0 = 1.;
    double h = 0.1;
    double t_end = 5.;
    
    Differential df;

    std::ofstream output("solution.csv");

    output << "t\teuler\trk4,exact\n";

    double t = t0;
    double y_euler = y0;
    double y_rk4 = y0;

    while (t <= t_end) {
        double exact = std::exp(-t);
        output << t << "\t" << y_euler << "\t" << y_rk4 << "\t" << exact << "\n" << std::endl;
        y_euler = df.euler_step(t, y_euler, h);
        y_rk4 = df.rk4_step(t, y_rk4, h);
        t += h;
    }
    output.close();
    std::cout << "done." << std::endl;
    return 0;
}