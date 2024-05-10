#include <iostream>
#include "MovingAverage.hpp"

int main(){
    MovingAverage<double> ma(3);

    std::cout << ma.step(1) << std::endl;
    std::cout << ma.step(10) << std::endl;
    std::cout << ma.step(3) << std::endl;
    std::cout << ma.step(5) << std::endl;
    std::cout << ma.step(7) << std::endl;

    ma.reset();

    std::cout << ma.get() << std::endl;

}