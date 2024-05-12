#include "Units.hpp"
#include "MovingAverage.hpp"
#include <vector>
#include<iostream>

class TrapezoidalMotionProfile{
    public:
    TrapezoidalMotionProfile(QSpeed maxVelocity, QAcceleration maxAcceleration) : maxVelocity(maxVelocity), maxAcceleration(maxAcceleration){
    }

    std::vector<QSpeed> generate(QLength distance){
        std::vector<QSpeed> profile((distance / maxVelocity).convert(10 * millisecond), maxVelocity);


        for(int i = 0; i < 1; i++){
            const std::size_t steps = (maxVelocity / maxAcceleration).convert(10 * millisecond);
            MovingAverage<QSpeed> average(steps);
            std::vector<QSpeed> newProfile;

            for(std::size_t _ = 0; _ < steps; _++){
                average.step(0_mps);
            }

            for(auto speed : profile){
                newProfile.push_back(average.step(speed));
            }

            for(std::size_t _ = 0; _ < steps; _++){
                newProfile.push_back(average.step(0_mps));
            }

            profile = newProfile;
        }

        return profile;

    }

    private:
    QSpeed maxVelocity;
    QAcceleration maxAcceleration;
};

class SCurveMotionProfile{
    public:
    SCurveMotionProfile(QSpeed maxVelocity, QAcceleration maxAcceleration, QJerk maxJerk) : 
        maxVelocity(maxVelocity), maxAcceleration(maxAcceleration), maxJerk(maxJerk){}

    std::vector<QSpeed> generate(QLength distance){
        std::vector<QSpeed> profile((distance / maxVelocity).convert(10 * millisecond), maxVelocity);

        for(int i = 0; i < 1; i++){
            const std::size_t steps = (maxVelocity / maxAcceleration).convert(10 * millisecond);
            MovingAverage<QSpeed> average(steps);
            std::vector<QSpeed> newProfile;

            for(std::size_t _ = 0; _ < steps; _++){
                average.step(0_mps);
            }

            for(auto speed : profile){
                newProfile.push_back(average.step(speed));
            }

            for(std::size_t _ = 0; _ < steps; _++){
                newProfile.push_back(average.step(0_mps));
            }

            profile = newProfile;
        }

        for(int i = 0; i < 1; i++){
            const std::size_t steps = (maxAcceleration / maxJerk).convert(10 * millisecond);
            MovingAverage<QSpeed> average(steps);
            std::vector<QSpeed> newProfile;

            for(std::size_t _ = 0; _ < steps; _++){
                average.step(0_mps);
            }

            for(auto speed : profile){
                newProfile.push_back(average.step(speed));
            }

            for(std::size_t _ = 0; _ < steps; _++){
                newProfile.push_back(average.step(0_mps));
            }

            profile = newProfile;
        }

        return profile;
    }

    private:
    QSpeed maxVelocity;
    QAcceleration maxAcceleration;
    QJerk maxJerk;
};