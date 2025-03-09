#pragma once

template<typename T>
class LinearFlux {
    T a;

public:
    LinearFlux(T a) : a(a) {}
    
    inline T operator() ([[maybe_unused]] T u) const {
        return a * u;
    }
};
