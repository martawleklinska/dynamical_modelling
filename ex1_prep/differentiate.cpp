#include "differentiate.hpp"
#include <iostream>
#include<fstream>


double Differential::func(double t, double y){
    return -y;
}

double Differential::euler_step(double t, double y, double h){
    return y + h * func(t, y);
}

double Differential::rk4_step(double t, double y, double h){
    double k1 = func(t, y);
    double k2 = func(t + h/2, y + k1 * h/2);
    double k3 = func(t + h/2, y + k2 * h/2);
    double k4 = func(t + h, y + h * k3);
    return y + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
}
