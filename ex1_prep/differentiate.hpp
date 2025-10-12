#pragma once
#include<iostream>
class Differential {
    public:
        double func(double t, double y);
        double euler_step(double t, double y, double h);
        double rk4_step(double t, double y, double h);
};