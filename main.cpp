#include <iostream>
#include "MovingAverage.hpp"
#include "MotionProfile.hpp"
#include "Units.hpp"

int main(){
    using namespace literals;

    SCurveMotionProfile profile(2_mps, 3_mps2, 8_mps3);
    std::vector<QSpeed> vels = profile.generate(6_m);

    for(auto vel : vels){
        std::cout << vel.convert(mps) << std::endl;
    }
}