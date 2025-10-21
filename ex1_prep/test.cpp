#include<iostream>
#include<cmath>

// Enum - pozwala wymienić liczby jako kolejne wyrazy

enum class colours {
    Blue,
    Red,
    Yellow,
    Green,
    ERROR=100,
};

// class Warrior {
//     private:

// };

int main(){
    colours col;
    col = colours::Blue;
    switch (col)
    {
    case colours::Blue:
        std::cout << "the color is blue" << std::endl;
        break;
    case colours::Yellow:
        std::cout << "the colour is yellow" << std::endl;
        break;
    case colours::Red:
        std::cout << "the colour is red" << std::endl;
        break;
    case colours::Green:
        std::cout << "the colour is green" << std::endl;
        break;
    // case colours::ERROR:
    //     std::cout << "the colour is error" << std::endl;
    //     break;
    default:
        break;
    }
    return 0;
}