#include "duffing.hpp"
#include <iostream>

int main() {
    Duffing duff(0.1, 1.0, -5.0, 0.0, 0.0);
    duff.solve();
    std::cout << "Symulacja zakończona.\n";
}
