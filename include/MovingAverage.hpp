#pragma once
#include <queue>

template<typename T>
class MovingAverage {
    public:
    MovingAverage(std::size_t size) : size(size) {}

    void reset(){
        while(!q.empty()){
            q.pop();
        }

        sum = T{};
    }

    T step(T val){
        if(q.size() == size){
            sum -= q.front();
            q.pop();
        }

        q.push(val);
        sum += val;

        return get();
    }

    T get() const{
        if(q.empty()){
            return T{};
        }

        return sum / q.size();
    }

    std::size_t getSize() const{
        return q.size();
    }

    std::size_t getCapacity() const{
        return size;
    }


    private:
    std::queue<T> q;
    std::size_t size;
    T sum{};
};