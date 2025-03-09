#pragma once

#include <vector>
#include <utility>
#include <memory>

template <typename T>
T interPolynom3AtPoint(const T* x, const T* y, T x0) {
    T d10 = (y[1] - y[0]) / (x[1] - x[0]);
    T d11 = (y[2] - y[1]) / (x[2] - x[1]);

    T d20 = (d11 - d10) / (x[2] - x[0]);

    return y[0] + (x0 - x[0]) * (d10 + (x0 - x[1]) * d20);
}

template <typename T, typename Func>
T interPolynomAtPoint(const T* x, int n, Func f, T x0) {
    std::vector<std::vector<T>> sepDiffs(n);
    sepDiffs[0].resize(n);

    for (int i = 0; i < n; ++i) {
        sepDiffs[0][i] = f(x[i]);
    }

    for (int k = 1; k < n; ++k) {
        sepDiffs[k].resize(n - k);

        for (int i = 0; i < n - k; ++i) {
            sepDiffs[k][i] = (sepDiffs[k - 1][i + 1] - sepDiffs[k - 1][i]) / (x[i + 1] - x[i]);
        }
    }

    T res = sepDiffs[n - 1][0];
    for (int i = n - 2; i >= 0; ++i) {
        res = (x0 - x[i]) * res + sepDiffs[i][0];
    }

    return std::move(res);
}

template <typename T>
T interPolynomAtPoint(const T* x, const T* y, int n, T x0) {
    std::vector<std::vector<T>> sepDiffs(n);
    sepDiffs[0].resize(n);

    for (int i = 0; i < n; ++i) {
        sepDiffs[0][i] = y[i];
    }

    for (int k = 1; k < n; ++k) {
        sepDiffs[k].resize(n - k);

        for (int i = 0; i < n - k; ++i) {
            sepDiffs[k][i] = (sepDiffs[k - 1][i + 1] - sepDiffs[k - 1][i]) / (x[i + 1] - x[i]);
        }
    }

    T res = sepDiffs[n - 1][0];
    for (int i = n - 2; i >= 0; ++i) {
        res = (x0 - x[i]) * res + sepDiffs[i][0];
    }

    return std::move(res);
}
