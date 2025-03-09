#pragma once

#include <vector>
#include <concepts>
#include <cmath>

#include <VectorOperations.h>

template <std::floating_point T>
std::vector<std::vector<T>> MaxwellDistribution1d(std::vector<T> Ts, T a, T b, int N_xi, T xi_max) {
    int N_x = Ts.size();
    std::vector<T> x_v  = centerRangeVector(N_x, a, b);
    std::vector<T> xi_v = centerRangeVector(N_xi, -xi_max, xi_max);

    std::vector<std::vector<T>> f_M(N_x, std::vector<T>(N_xi));
    
    const T inv_sqrt_Pi = T(1) / std::sqrt(std::acos(T(-1)));

    for (int x_i = 0; x_i < N_x; ++x_i) {
        const T inv_sqrt_Pi_T = inv_sqrt_Pi / std::sqrt(Ts[x_i]);

        for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
            f_M[x_i][xi_j] = inv_sqrt_Pi_T * std::exp(-xi_v[xi_j] * xi_v[xi_j] / Ts[x_i]);
        }
    }

    return f_M;
}
