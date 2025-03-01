#pragma once

#include <memory>

/*
// Rusanov monotone flux struct
template <typename T>
struct RusanovFlux {
    static inline T flux(T fl, T fr, T ul, T ur) {
        return T(0.5f) * (fl + fr - std::max(std::fabs(ur), std::fabs(ul)) * (ur - ul));
    }
};
*/

/*
// Rusanov monotone flux class (singleton)
template <typename T>
class RusanovFlux {
private:
    static std::unique_ptr<RusanovFlux<T>> instance;
    
    //template <typename C>
    // friend std::unique_ptr<RusanovFlux<T>> std::make_unique<RusanovFlux<T>>(); // Doesn't work
    RusanovFlux() noexcept {}

public:
    static RusanovFlux<T>* const getInstance();

    inline T operator() (T fl, T fr, T ul, T ur) const {
        return T(0.5f) * (fl + fr - std::max(std::fabs(ur), std::fabs(ul)) * (ur - ul));
    }
};

template <typename T>
RusanovFlux<T>* const RusanovFlux<T>::getInstance() {
    if (!instance) {
        // instance = std::make_unique<RusanovFlux<T>>(); // Best way to do
        instance = std::unique_ptr<RusanovFlux<T>>(new RusanovFlux());
    }

    return instance.get();
}
*/

// Rusanov monotone flux
template <typename T>
inline T rusanovFlux(T fl, T fr, T ul, T ur) {
    return T(0.5f) * (fl + fr - std::max(std::fabs(ur), std::fabs(ul)) * (ur - ul));
}

// Lax-Friedrichs monotone flux for uniform mesh
template <typename T>
class UniformSimpleLaxFriedrichsFlux {
    T _deltaX;
    T _deltaT;
    T _halfS;

public:
    UniformSimpleLaxFriedrichsFlux(T deltaX, T deltaT)
        : _deltaX(deltaX), _deltaT(deltaT), _halfS(T(0.5f)* deltaX / deltaT)
    {}

    inline T operator()(T fl, T fr, T ul, T ur) const {
        return T(0.5f) * (fl + fr) - _halfS * (ur - ul);
    }
};

// Accurate flux for constant wind direction
template <typename T>
class AccurateFlux {
    bool isRightwind; // is true, if wind is directed to right

public:
    AccurateFlux(bool isRightwind) : isRightwind(isRightwind) {}

    inline T operator()(T fl, T fr, [[maybe_unused]] T ul, [[maybe_unused]] T ur) const {
        return isRightwind * fl + !isRightwind * fr;
    }
};
