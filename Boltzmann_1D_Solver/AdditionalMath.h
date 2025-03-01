#pragma once

#include <cmath>

template <typename T>
inline constexpr int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
inline constexpr T minmod(T a, T b) {
    return T(0.5f) * (sign(a) + sign(b)) * std::fmin(std::fabs(a), std::fabs(b));
}

template <typename T>
inline constexpr T sqr(T x) { // Hope inline copies x for one evaluation, otherwise not optimised
    return x * x;
}
