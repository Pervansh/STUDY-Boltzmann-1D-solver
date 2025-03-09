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