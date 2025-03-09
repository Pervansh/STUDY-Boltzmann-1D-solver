#pragma once

#include <concepts>

/*
    концепт для RK: rkMethod::step(...) -> Type
*/
/*
template <
    typename S,
    std::floating_point Type,
    typename Y,
    typename Func
>
*/

/*
    концепт для Func: Func(Type u, Type t) -> Type
*/

template <std::floating_point T, typename Y, typename F>
using RkStepY = Y(*)(F f, Y y, T tau);

// Type - time type, Y - phase coordinate type
struct SSPRK3 {
    /*
    static unsigned int svtages;
    static Type a[];
    static Type b[];
    */

    // RK step for f w/ signature f(Type, Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepTY(Func f, T t, Y y, T tau) {
        Y k1 = f(t, y);
        Y k2 = f(t + tau, y + tau * k1);
        Y k3 = f(t + T(0.5f) * tau, y + T(0.25f) * tau * (k1 + k2));

        return y + tau / T(6.f) * (k1 + k2 + T(4.f) * k3);
    }

    // RK step for f w/ signature f(Y, Type)
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

// Type - time type, Y - phase coordinate type
struct HeunsMethodRK {
    /*
    static unsigned int stages;
    static Type a[];
    static Type b[];
    */
    // RK step for f w/ signature f(Type, Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepTY(Func f, T t, Y y, T tau) {
        Y k1 = f(t, y);
        Y k2 = f(t + tau, y + tau * k1);

        return y + tau / T(2.f) * (k1 + k2);
    }

    // RK step for f w/ signature f(Y, Type)
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


// Type - time type, Y - phase coordinate type
struct ExplicitEulerRK {
    /*
    static unsigned int stages;
    static Type a[];
    static Type b[];
    */
    // RK step for f w/ signature f(Type, Y)
    template<typename T, typename Y, typename Func>
    static inline Y stepTY(Func f, T t, Y y, T tau) {
        return y + tau * f(t, y);
    }

    // RK step for f w/ signature f(Y, Type)
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
template <typename Type>
unsigned int SSPRK3<Type>::stages = 3;

template <typename Type>
Type SSPRK3<Type>::a[] = {Type(0), Type(1), Type(0.5)};

template <typename Type>
Type SSPRK3<Type>::b[] = { 
    {Type(0), Type(1), Type(0.5)},
    {Type(1), Type(0), Type(0)}
    {Type(0.25), Type(0.25), Type(0)}
};
*/
/*
template<typename Type, typename Func>
vector<Type> rungeKuttaTemplateStep(const RungeKuttaParams<Type>& rkParams, Type tau, int n, Func f, Type t, const vector<Type>& y) {
    int m = rkParams.stages;
    vector<vector<Type>> k(m);

    auto& a = rkParams.a;
    auto& b = rkParams.b;
    auto& sigma = rkParams.sigma;

    k[0] = f(t, y);

    vector<Type> kSum = vector<Type>(n, 0);
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



