#pragma once

#include <vector>
#include <concepts>

template<typename F, typename T>
concept ConstIntegrant = std::floating_point<T> && requires (const F& func, T val) {
    { func(val) } -> std::same_as<T>;
};

template<typename F, typename T>
concept ApproxInt = std::floating_point<T> && requires (const F & func, const std::vector<T> ys, T dx) {
    { func(ys, dx) } -> std::same_as<T>;
};

template<typename T>
using ApproxIntFunc = T(*)(const std::vector<T> ys, T dx);


// Declarations

template<std::floating_point T>
T Cell1stOrderInt(const std::vector<T> ys, T dx);

template<std::floating_point T>
T simpsonInt(const std::vector<T> ys, T dx);

// Definitions

template<std::floating_point T>
T Cell1stOrderInt(const std::vector<T> ys, T dx) {
    T res{};
    int N = ys.size();

    for (int i = 0; i < N; ++i) {
        res += dx * ys[i];
    }

    return res;
}

template<std::floating_point T>
T simpsonInt(const std::vector<T> ys, T dx) {
    T res{};
    int N = ys.size();

    if (N <= 0) {
        return T(0.);
    }

    if (N == 1) {
        return dx * ys[0];
    }

    const T dxThird = dx / T(3.);

    for (int i = 0; i + 2 < N; i += 2) {
        res += dxThird * (ys[i] + 4 * ys[i + 1] + ys[i + 2]);
    }

    if (N % 2 == 0) {
        res += dx * (ys[N - 2] + ys[N - 1]) / T(2.);
    }

    return res;
}