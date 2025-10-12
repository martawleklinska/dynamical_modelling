#include "oscillator.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

std::pair<double, double> Harmonic::func(double x, double y){
    double dx = y;
    double dy = -x;
    return {dx, dy};
}


int main() {
    double xmin = -2.0, xmax = 2.0;
    double ymin = -2.0, ymax = 2.0;
    int nx = 20, ny = 20;

    Harmonic harm;

    std::ofstream out("vector_field.csv");
    out << "x,y,dx,dy\n";

    for (int i = 0; i < nx; ++i) {
        double x = xmin + i * (xmax - xmin) / (nx - 1);
        for (int j = 0; j < ny; ++j) {
            double y = ymin + j * (ymax - ymin) / (ny - 1);
            auto [dx, dy] = harm.func(x, y);
            out << x << "," << y << "," << dx << "," << dy << "\n";
        }
    }
    out.close();
    std::cout << "Saved vector field to vector_field.csv\n";
}
