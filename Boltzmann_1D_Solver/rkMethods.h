#pragma once

#include <type_traits>

/*
    концепт для RK: rkMethod::step(...) -> T
*/

/*
    концепт для Func: Func(T u, T t) -> T
*/

// T - time type, Y - phase coordinate type
struct SSPRK3 {
    /*
    static unsigned int svtages;
    static T a[];
    static T b[];
    */

    // RK step for f w/ signature f(T, Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepTY(Func f, T t, Y y, T tau) {
        Y k1 = f(t, y);
        Y k2 = f(t + tau, y + tau * k1);
        Y k3 = f(t + T(0.5f) * tau, y + T(0.25f) * tau * (k1 + k2));

        return y + tau / T(6.f) * (k1 + k2 + T(4.f) * k3);
    }

    // RK step for f w/ signature f(Y, T)
    template<typename T, typename Y, typename Func>
    static inline Y stepYT(Func f, Y y, T t, T tau) {
        Y k1 = f(y, t);
        Y k2 = f(y + tau * k1, t + tau);
        Y k3 = f(y + T(0.25f) * tau * (k1 + k2), t + T(0.5f) * tau);

        return y + tau / T(6.f) * (k1 + k2 + T(4.f) * k3);
    }

    // RK step for f w/ signature f(Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepY(Func f, Y y, T tau) {
        Y k1 = f(y);
        Y k2 = f(y + tau * k1);
        Y k3 = f(y + T(0.25f) * tau * (k1 + k2));

        return y + tau / T(6.f) * (k1 + k2 + T(4.f) * k3);
    }
};

// T - time type, Y - phase coordinate type
struct HeunsMethodRK {
    /*
    static unsigned int stages;
    static T a[];
    static T b[];
    */
    // RK step for f w/ signature f(T, Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepTY(Func f, T t, Y y, T tau) {
        Y k1 = f(t, y);
        Y k2 = f(t + tau, y + tau * k1);

        return y + tau / T(2.f) * (k1 + k2);
    }

    // RK step for f w/ signature f(Y, T)
    template<typename T, typename Y, typename Func>
    static inline Y stepYT(Func f, Y y, T t, T tau) {
        Y k1 = f(y, t);
        Y k2 = f(y + tau * k1, t + tau);

        return y + tau / T(2.f) * (k1 + k2);
    }

    // RK step for f w/ signature f(Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepY(Func f, Y y, T tau) {
        Y k1 = f(y);
        Y k2 = f(y + tau * k1);

        return y + tau / T(2.f) * (k1 + k2);
    }
};


// T - time type, Y - phase coordinate type
struct ExplicitEulerRK {
    /*
    static unsigned int stages;
    static T a[];
    static T b[];
    */
    // RK step for f w/ signature f(T, Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepTY(Func f, T t, Y y, T tau) {
        return y + tau * f(t, y);
    }

    // RK step for f w/ signature f(Y, T)
    template<typename T, typename Y, typename Func>
    static inline Y stepYT(Func f, Y y, T t, T tau) {
        return y + tau * f(y, t);
    }

    // RK step for f w/ signature f(Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepY(Func f, Y y, T tau) {
        return y + tau * f(y);
    }
};

/*
template <typename T>
unsigned int SSPRK3<T>::stages = 3;

template <typename T>
T SSPRK3<T>::a[] = {T(0), T(1), T(0.5)};

template <typename T>
T SSPRK3<T>::b[] = { 
    {T(0), T(1), T(0.5)},
    {T(1), T(0), T(0)}
    {T(0.25), T(0.25), T(0)}
};
*/
/*
template<typename T, typename Func>
vector<T> rungeKuttaTemplateStep(const RungeKuttaParams<T>& rkParams, T tau, int n, Func f, T t, const vector<T>& y) {
    int m = rkParams.stages;
    vector<vector<T>> k(m);

    auto& a = rkParams.a;
    auto& b = rkParams.b;
    auto& sigma = rkParams.sigma;

    k[0] = f(t, y);

    vector<T> kSum = vector<T>(n, 0);
    for (int i = 1; i < m; i++) {
        kSum.assign(n, 0);
        for (int j = 0; j < i; ++j) {
            kSum += b[i][j] * k[j];
        }

        k[i] = f(t + a[i] * tau, y + tau * kSum);
    }

    kSum.assign(n, 0);
    for (int j = 0; j < m; ++j) {
        kSum += sigma[j] * k[j];
    }

    auto nextY = y + tau * kSum;

    return nextY;
}
*/



