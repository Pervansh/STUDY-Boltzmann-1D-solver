#pragma once

#include <cmath>
#include <concepts>

constexpr double D_PI = 3.14159265358979323846;

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

//*
// \returns pos_val when cond_val > 0, neg_val when cond_val < 0, (pos_val + neg_val)/2 otherwise
//
template <std::floating_point T, typename U = T>
requires requires(U u) { sign(u); }
inline constexpr T signCond(T pos_val, T neg_val, U cond_val) {
    T cond_sgn(sign(cond_val));
    return T(0.5f) * ((1 + cond_sgn) * pos_val + (1 - cond_sgn) * neg_val);
}
