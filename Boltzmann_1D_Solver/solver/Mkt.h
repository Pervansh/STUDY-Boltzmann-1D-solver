#pragma once

#include <vector>
#include <concepts>
#include <cmath>

#include <VectorOperations.h>

template<std::floating_point T>
using ViscosityFunc = T(*)(T Temperature);

template<std::floating_point T>
inline T ballMoleculeViscosityRule(T Temp) {
    return std::sqrt(Temp);
}

template<std::floating_point T>
inline T maxwellMoleculeViscosityRule(T Temp) {
    return Temp;
}

template<std::floating_point T>
std::vector<T> normalDistribution1d(T sigma, int N, T x_max) {
    std::vector<T> xi_v = centerRangeVector(N, -x_max, x_max);

    std::vector<T> f(N);

    const T inv_sqrt_Pi = T(1) / std::sqrt(2 * std::acos(T(-1)));

    const T doubled_sigma_sqr = 2 * sigma * sigma;

#pragma omp parallel for
    for (int xi_j = 0; xi_j < N; ++xi_j) {
        f[xi_j] = inv_sqrt_Pi * std::exp(-xi_v[xi_j] * xi_v[xi_j] / doubled_sigma_sqr) / sigma; 
    }

    return f;
}

template <std::floating_point T>
std::vector<std::vector<T>> MaxwellDistribution1d(std::vector<T> Ts, T a, T b, int N_xi, T xi_max) {
    int N_x = Ts.size();
    std::vector<T> x_v  = centerRangeVector(N_x, a, b);
    std::vector<T> xi_v = centerRangeVector(N_xi, -xi_max, xi_max);

    std::vector<std::vector<T>> f_M(N_x, std::vector<T>(N_xi));
    
    const T inv_sqrt_Pi = T(1) / std::sqrt(std::acos(T(-1)));

#pragma omp parallel for
    for (int x_i = 0; x_i < N_x; ++x_i) {
        const T inv_sqrt_Pi_T = inv_sqrt_Pi / std::sqrt(Ts[x_i]);
        auto& f_M_cur_xi = f_M[x_i];

        for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
            f_M_cur_xi[xi_j] = inv_sqrt_Pi_T * std::exp(-xi_v[xi_j] * xi_v[xi_j] / Ts[x_i]);
        }
    }

    return f_M;
}

template <std::floating_point T>
std::vector<T> MaxwellDistribution1d(T Temp, int N_xi, T xi_max) {
    std::vector<T> xi_v = centerRangeVector(N_xi, -xi_max, xi_max);

    std::vector<T> f_M(N_xi);

    const T inv_sqrt_Pi = T(1) / std::sqrt(std::acos(T(-1)));

    const T inv_sqrt_Pi_T = inv_sqrt_Pi / std::sqrt(Temp);

#pragma omp parallel for
    for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
        f_M[xi_j] = inv_sqrt_Pi_T * std::exp(-xi_v[xi_j] * xi_v[xi_j] / Temp);
    }

    return f_M;
}

template <std::floating_point T>
std::vector<T> unnormalizedMaxwellDistribution1d(T Temp, int N_xi, T xi_max) {
    std::vector<T> xi_v = centerRangeVector(N_xi, -xi_max, xi_max);

    std::vector<T> f_M(N_xi);

    const T inv_sqrt_Pi = T(1) / std::sqrt(std::acos(T(-1)));

#pragma omp parallel for
    for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
        f_M[xi_j] = inv_sqrt_Pi * std::exp(-xi_v[xi_j] * xi_v[xi_j] / Temp);
    }

    return f_M;
}
