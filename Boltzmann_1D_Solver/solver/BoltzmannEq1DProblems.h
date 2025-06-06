#pragma once

#include <vector>
#include <assert.h>
#include <concepts>
#include <cmath>

#include "BoltzmannEq1DSolvers.h"
#include "Mkt.h"

template<std::floating_point T>
Bgk1dProblemData<T> temperatureRiemannProblem(T T_1, T T_2, int N_x, T x_max, int N_xi, T xi_max, T dt, T t_end, T delta) {
    assert(N_x % 2 == 0);

    Bgk1dProblemData<T> data;
    data.N_x = N_x;
    data.N_xi = N_xi;
    data.dt = dt;
    data.t_end = t_end;
    data.xi_max = xi_max;
    data.delta = delta;

    data.a = -x_max;
    data.b = x_max;

    data.Ts = std::vector<T>(N_x);

    for (int x_i = 0; x_i < N_x / 2; ++x_i) {
        data.Ts[x_i] = T_1;
    }

    for (int x_i = N_x / 2; x_i < N_x; ++x_i) {
        data.Ts[x_i] = T_2;
    }

    data.f_1 = MaxwellDistribution1d<T>(data.Ts, data.a, data.b, N_xi, data.xi_max);

    return data;
}

template<std::floating_point T>
Bgk1dProblemData<T> densityRiemannProblem(
    T Temp, T n_1, T n_2,
    int N_x, T x_max, int N_xi, T xi_max,
    T dt, T t_end,
    T delta,
    T discont_point = T(0)
) {
    assert(N_x % 2 == 0);

    Bgk1dProblemData<T> data;
    data.N_x = N_x;
    data.N_xi = N_xi;
    data.dt = dt;
    data.t_end = t_end;
    data.xi_max = xi_max;
    data.delta = delta;

    data.a = -x_max;
    data.b = x_max;

    T dx = 2 * x_max / N_x;

    int discont_x_ind = (int)std::round((discont_point - data.a) / dx);

    data.Ts = std::vector<T>(N_x);

    for (int x_i = 0; x_i < N_x; ++x_i) {
        data.Ts[x_i] = Temp;
    }

    data.f_1 = MaxwellDistribution1d<T>(data.Ts, data.a, data.b, N_xi, data.xi_max);

    for (int x_i = 0; x_i < discont_x_ind; ++x_i) {
        data.f_1[x_i] *= n_1;
    }

    for (int x_i = discont_x_ind; x_i < N_x; ++x_i) {
        data.f_1[x_i] *= n_2;
    }

    return data;
}

template<std::floating_point T>
Bgk1dProblemData<T> densityRiemannProblem(
    T Temp, T n_1, T n_2,
    int N_x, T  a, T b, int N_xi, T xi_max,
    T dt, T t_end,
    T delta,
    T discont_point
) {
    assert(N_x % 2 == 0);

    Bgk1dProblemData<T> data;
    data.N_x = N_x;
    data.N_xi = N_xi;
    data.dt = dt;
    data.t_end = t_end;
    data.xi_max = xi_max;
    data.delta = delta;

    data.a = a;
    data.b = b;

    T dx = (b - a) / N_x;

    int discont_x_ind = (int)std::round((discont_point - data.a) / dx);

    data.Ts = std::vector<T>(N_x);

    for (int x_i = 0; x_i < N_x; ++x_i) {
        data.Ts[x_i] = Temp;
    }

    data.f_1 = MaxwellDistribution1d<T>(data.Ts, data.a, data.b, N_xi, data.xi_max);

    for (int x_i = 0; x_i < discont_x_ind; ++x_i) {
        data.f_1[x_i] *= n_1;
    }

    for (int x_i = discont_x_ind; x_i < N_x; ++x_i) {
        data.f_1[x_i] *= n_2;
    }

    return data;
}

template<std::floating_point T>
Bgk1dProblemData<T> densityRiemannProblem(
    T Temp, T n_1, T n_2,
    int N_x, T  a, T b, int N_xi, T xi_max,
    T dt, T t_end,
    T delta
) {
    return densityRiemannProblem<T>(Temp, n_1, n_2, N_x, a, b, N_xi, xi_max, dt, t_end, delta, a + b / 2);
}

template<std::floating_point T>
Bgk1dProblemData<T> densityTemperatureRiemannProblem(
    T Temp_1, T n_1,
    T Temp_2, T n_2,
    int N_x, T  a, T b, int N_xi, T xi_max,
    T dt, T t_end,
    T delta,
    T discont_point
) {
    assert(N_x % 2 == 0);

    Bgk1dProblemData<T> data;
    data.N_x = N_x;
    data.N_xi = N_xi;
    data.dt = dt;
    data.t_end = t_end;
    data.xi_max = xi_max;
    data.delta = delta;

    data.a = a;
    data.b = b;

    T dx = (b - a) / N_x;

    int discont_x_ind = (int)std::round((discont_point - data.a) / dx);

    data.Ts = std::vector<T>(N_x);

    for (int x_i = 0; x_i < discont_x_ind; ++x_i) {
        data.Ts[x_i] = Temp_1;
    }

    for (int x_i = discont_x_ind; x_i < N_x; ++x_i) {
        data.Ts[x_i] = Temp_2;
    }

    data.f_1 = MaxwellDistribution1d<T>(data.Ts, data.a, data.b, N_xi, data.xi_max);

    for (int x_i = 0; x_i < discont_x_ind; ++x_i) {
        data.f_1[x_i] *= n_1;
    }

    for (int x_i = discont_x_ind; x_i < N_x; ++x_i) {
        data.f_1[x_i] *= n_2;
    }

    return data;
}

template<std::floating_point T>
Bgk1dProblemData<T> densityTemperatureRiemannProblem(
    T Temp_1, T n_1,
    T Temp_2, T n_2,
    int N_x, T  a, T b, int N_xi, T xi_max,
    T dt, T t_end,
    T delta
) {
    return densityTemperatureRiemannProblem<T>(Temp_1, n_1, Temp_2, n_2, N_x, a, b, N_xi, xi_max, dt, t_end, delta, a + b / 2);
}

/*
    Returns problem data for emitting wall problem
    Region: 0 < x < x_max, 0 < t < t_end
    I.C.: n = n_0, T = T_0
    B.C.: x = 0, xi > 0: h(x, xi) = h_M(xi, n_wall, T_wall), g(x, xi) = g_M(xi, n_wall, T_wall)
*/
template<std::floating_point T>
Bgk1dProblemData<T> emittingWallProblem(
    T Temp_0, T n_0,
    T Temp_wall, T n_wall,
    int N_x, T x_max, int N_xi, T xi_max,
    T dt, T t_end,
    T delta
) {
    Bgk1dProblemData<T> data;
    data.N_x = N_x;
    data.N_xi = N_xi;
    data.dt = dt;
    data.t_end = t_end;
    data.xi_max = xi_max;
    data.delta = delta;

    data.a = T(0.);
    data.b = x_max;

    T dx = x_max / N_x;

    // int discont_x_ind = (int)std::round((discont_point - data.a) / dx);

    data.Ts = std::vector<T>(N_x);

    for (int x_i = 0; x_i < N_x; ++x_i) {
        data.Ts[x_i] = Temp_0;
    }

    data.f_1 = MaxwellDistribution1d<T>(data.Ts, data.a, data.b, N_xi, data.xi_max);

    std::vector<T> f_wall = MaxwellDistribution1d(Temp_wall, N_xi, xi_max);
    data.f_1.front() = std::move(f_wall);

    data.f_1.front() *= n_wall;

    for (int x_i = 1; x_i < N_x; ++x_i) {
        data.f_1[x_i] *= n_0;
    }

    return data;
}

/*
    Returns problem data for emitting wall problem
    Region: 0 < x < x_max, 0 < t < t_end
    I.C.: n ≡ n_0, T ≡ T_0
    B.C.: x = 0, xi > 0: h(x, xi) = h_M(xi, n_wall, T_wall), g(x, xi) = g_M(xi, n_wall, T_wall)
*/
template<std::floating_point T>
Bgk1dProblemData<T> evaporatingWallProblem(
    T Temp_0, T n_0,
    T Temp_wall, T Q_wall,
    int N_x, T x_max, int N_xi, T xi_max,
    T dt, T t_end,
    T delta
) {
    Bgk1dProblemData<T> data;
    data.N_x = N_x;
    data.N_xi = N_xi;
    data.dt = dt;
    data.t_end = t_end;
    data.xi_max = xi_max;
    data.delta = delta;

    data.a = T(0.);
    data.b = x_max;

    T dx = x_max / N_x;

    // int discont_x_ind = (int)std::round((discont_point - data.a) / dx);

    data.Ts = std::vector<T>(N_x);

    for (int x_i = 0; x_i < N_x; ++x_i) {
        data.Ts[x_i] = Temp_0;
    }

    data.f_1 = MaxwellDistribution1d<T>(data.Ts, data.a, data.b, N_xi, data.xi_max);

    std::vector<T> f_wall = MaxwellDistribution1d(Temp_wall, N_xi, xi_max); // original
    // std::vector<T> f_wall = unnormalizedMaxwellDistribution1d(Temp_wall, N_xi, xi_max); // bruteforce

    T coef = std::exp(-Q_wall * (T(1.) / Temp_wall - T(1.))) / Temp_wall;

    // f_wall *= std::exp(-Q_wall * (T(1.) / Temp_wall - T(1.))) / std::sqrt(Temp_wall); // correct?
    // f_wall *= std::exp(-Q_wall * (T(1.) / Temp_wall - T(1.))) / Temp_wall; // incorrect
    // f_wall *= 74. / std::sqrt(D_PI); // bruteforce possible correction
    f_wall *= std::exp(-Q_wall * (T(1.) / Temp_wall - T(1.))) / Temp_wall;

    auto xi_v = centerRangeVector(N_xi, -xi_max, xi_max);

    for (int i = 0; i < N_xi; ++i) {
        if (xi_v[i] <= 0.) {
            f_wall[i] = 0.;
        }
    }

    data.f_1.front() = std::move(f_wall);

    for (int x_i = 1; x_i < N_x; ++x_i) {
        data.f_1[x_i] *= n_0;
    }

    return data;
}
