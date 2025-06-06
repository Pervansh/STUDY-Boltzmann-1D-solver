#pragma once

#include <vector>
#include <iostream>
#include <utility>
#include <concepts>
#include <cmath>

#include <assert.h>

#include "AdditionalMath.h"
#include "Norma.h"
#include "rkMethods.h"
#include "VectorOperations.h"
#include "Integration.h"
#include "MacroparametersFileInterface.h"

template<typename Func, typename T>
concept OutputRuleType = std::floating_point<T> && requires (Func rule, T t) { { rule(t) } -> std::same_as<bool>; };

template<std::floating_point T>
T knToDelta(T Kn) {
    return 8. / (5. * D_PI * Kn);
}

template <std::floating_point T>
struct Destribution1dState{
    std::vector<std::vector<T>> h; // h[i][j] = h_ij
    std::vector<std::vector<T>> g; // g[i][j] = g_ij

    Destribution1dState() = default;

    // Perfect Forwarding
    template<
        std::convertible_to< std::vector<std::vector<T>> > H,
        std::convertible_to< std::vector<std::vector<T>> > G
    >
    Destribution1dState(H&& h_, G&& g_) : h(std::forward<H>(h_)), g(std::forward<G>(g_)) {}

    Destribution1dState(int N_x, int N_xi) :
        h(N_x, std::vector<T>(N_xi)),
        g(N_x, std::vector<T>(N_xi))
    {}
};

template <std::floating_point T>
inline Destribution1dState<T> operator+(
    const Destribution1dState<T>& s1,
    const Destribution1dState<T>& s2
) {
    return Destribution1dState<T>(s1.h + s2.h, s1.g + s2.g);
}

template <std::floating_point T>
inline Destribution1dState<T> operator*(T k, const Destribution1dState<T>& s) {
    return Destribution1dState<T>(k * s.h, k * s.g);
}

template <std::floating_point T>
struct Bgk1dProblemData {
    int N_x;
    int N_xi;
    T a;
    T b;
    T xi_max;
    T dt;
    T t_end;
    std::vector<T> Ts; // Ts[i] = T_i - avg. temperature in i-th cell
    std::vector<std::vector<T>> f_1; // f_1(x_1, xi_1)
    T delta;
};

template <std::floating_point Type>
struct Macroparameters1d {
    std::vector<Type> n;
    std::vector<Type> u_1;
    std::vector<Type> T;
    std::vector<Type> q_1;

    Macroparameters1d(int N_x) {
        n.resize(N_x);
        u_1.resize(N_x);
        T.resize(N_x);
        q_1.resize(N_x);
    }
};

//*
// \param xi_v should be monotonous
//
template <std::floating_point Type>
Macroparameters1d<Type> bgk1dMacroparameters(
    const Destribution1dState<Type>& state,
    const Type dxi,
    const std::vector<Type>& xi_v,
    const ApproxIntFunc<Type>& apprInt
) {
    assert(state.h.size() == state.g.size());
    int N_x = state.h.size();

    Macroparameters1d<Type> params(N_x);

    for (int x_i = 0; x_i < N_x; ++x_i) {
        assert(state.h[x_i].size() == state.g[x_i].size());

        int N_xi = state.h[x_i].size();
        assert(N_xi == xi_v.size());

        Type& n = params.n[x_i];
        Type& u_1 = params.u_1[x_i];
        Type& T = params.T[x_i];
        Type& q_1 = params.q_1[x_i];

        std::vector<Type> u_1_int_f(N_xi);
        std::vector<Type> T_int_f(N_xi);
        std::vector<Type> q_1_int_f(N_xi);

        // Integrant calculations for u_1 & Type
        for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
            u_1_int_f[xi_j] = xi_v[xi_j] * state.h[x_i][xi_j];
            T_int_f[xi_j] = xi_v[xi_j] * xi_v[xi_j] * state.h[x_i][xi_j] + state.g[x_i][xi_j];
        }

        // Macroparameters calculations for n, u_1 & Type
        n = apprInt(state.h[x_i], dxi);
        u_1 = apprInt(u_1_int_f, dxi) / n;
        Type T_int = apprInt(T_int_f, dxi);
        T = 2. / 3. * (T_int / n - u_1 * u_1);

        // Integrant calculations for q_1
        for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
            Type v = xi_v[xi_j] - u_1;
            q_1_int_f[xi_j] = v * (v * v * state.h[x_i][xi_j] + state.g[x_i][xi_j]);
        }

        // Macroparameters calculations for q_1 & c_k
        q_1 = 0.5 * apprInt(q_1_int_f, dxi);
    }

    return params;
}

template <typename RkMethod, std::floating_point T>
inline void bgk1dMethod(
    const Bgk1dProblemData<T>& data,
    const ApproxIntFunc<T>& apprInt,
    Full1dStateOutput& output,
    const ViscosityFunc<T>& mu = ballMoleculeViscosityRule<T>
) {
    bgk1dMethod<RkMethod>(data, apprInt, output, [](T) {return false; }, mu);
}

template <typename RkMethod, std::floating_point T, OutputRuleType<T> OutputRule>
void bgk1dMethod(
    const Bgk1dProblemData<T>& data,
    const ApproxIntFunc<T>& apprInt,
    Full1dStateOutput& output,
    const OutputRule& outputRule,
    const ViscosityFunc<T>& mu = ballMoleculeViscosityRule<T>
) {
    int N_x = data.N_x;
    int N_xi = data.N_xi;

    assert(data.f_1.size() == N_x);

    std::vector<std::vector<T>> g(N_x);
    
    for (int x_i = 0; x_i < N_x; ++x_i) {
        assert(data.f_1[x_i].size() == N_xi);
        g[x_i].resize(N_xi);

        for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
            g[x_i][xi_j] = data.Ts[x_i] * data.f_1[x_i][xi_j];
        }
    }

    Destribution1dState<T> state;
    state.h = data.f_1;
    state.g = std::move(g);

    T dx = T{ (data.b - data.a) } / N_x;
    T dxi = T{ 2 * data.xi_max } / N_xi;

    std::vector<T> x_v(N_x);
    std::vector<T> xi_v(N_xi);

    for (int x_i = 0; x_i < N_x; ++x_i) {
        x_v[x_i] = data.a + (x_i + T(0.5f)) * dx;
    }

    for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
        xi_v[xi_j] = -data.xi_max + (xi_j + T(0.5f)) * dxi;
    }

    // Consts
    const T inv_sqrt_Pi = T(1.) / std::sqrt(std::acos(T(-1.)));

    // Boundry conditions arrays
    std::vector<T> h_l = state.h[0];
    std::vector<T> h_r = state.h[N_x - 1];
    std::vector<T> g_l = state.g[0];
    std::vector<T> g_r = state.g[N_x - 1];

    // Universal right side lambdas & calculation arrays
    std::vector<T> fL_j_v(N_x + 1);
    std::vector<T> fR_j_v(N_x + 1);
    std::vector<T> G_j_v(N_x + 1);

    std::vector<std::vector<T>> J_h(N_x, std::vector<T>(N_xi));
    std::vector<std::vector<T>> J_g(N_x, std::vector<T>(N_xi));

    auto xisRightSide(
        [&N_x, &fL_j_v, &fR_j_v, &G_j_v, &xi_v, &dx] (
            const std::vector<std::vector<T>>& f_v,
            const std::vector<std::vector<T>>& J_v,
            int xi_j, 
            std::vector<std::vector<T>>& diffs,
            const T f_l,
            const T f_r
        ) {
            // Domain 1st order boundry conditions and treatment (equilibrium medium / vacuum)

            T halfDelta0 = T(0.5f) * minmod(f_v[0][xi_j] - f_l, f_v[1][xi_j] - f_v[0][xi_j]); // 1st order condition
            fL_j_v[1] = f_v[0][xi_j] + halfDelta0; // Left value for right cell boundry
            fR_j_v[0] = f_v[0][xi_j] - halfDelta0; // Right value for left cell boundry
            // fR_j_v[N_x] = fR_j_v[0]; // periodic condition
            fR_j_v[N_x] = f_r; // 1st order condition

            T halfDeltaN_x = T(0.5f)* minmod(f_v[N_x - 1][xi_j] - f_v[N_x - 2][xi_j], f_r - f_v[N_x - 1][xi_j]); // 1st order condition
            fL_j_v[N_x] = f_v[N_x - 1][xi_j] + halfDeltaN_x; // Left value for right cell boundry
            fR_j_v[N_x - 1] = f_v[N_x - 1][xi_j] - halfDeltaN_x; // Right value for left cell boundry
            // fL_j_v[0] = fL_j_v[N_x]; // periodic condition
            fL_j_v[0] = f_l; // 1st order condition

            // Calculations inner cell boundries (iterating over cells)
            for (int i = 1; i < N_x - 1; ++i) {
                T halfDelta = T(0.5f) * minmod(f_v[i][xi_j] - f_v[i - 1][xi_j], f_v[i + 1][xi_j] - f_v[i][xi_j]);
                fL_j_v[i + 1] = f_v[i][xi_j] + halfDelta; // Left value for right cell boundry
                fR_j_v[i] = f_v[i][xi_j] - halfDelta; // Right value for left cell boundry
            }

            // Monotone fluxes calculation for all cell's boundries
            T cur_xi = xi_v[xi_j];
            for (int i = 0; i <= N_x; ++i) {
                G_j_v[i] = cur_xi * signCond(fL_j_v[i], fR_j_v[i], cur_xi);
            }

            for (int i = 0; i < N_x; ++i) {
                diffs[i][xi_j] = (G_j_v[i] - G_j_v[i + 1]) / dx + J_v[i][xi_j];
            }
        }
    );

    auto fullRightSide(
        [
            &N_x, &N_xi, &dx, dxi, &xi_v, &h_l, &h_r, &g_l, &g_r, &J_h, &J_g, 
            &apprInt, &mu, &data, &xisRightSide, inv_sqrt_Pi
        ] (
        const Destribution1dState<T>& state
    ) {
        // Should be passed through lambda argument (no dynamic init here)
        Destribution1dState<T> diff_state(N_x, N_xi);

        Macroparameters1d<T> params = bgk1dMacroparameters(state, dxi, xi_v, apprInt);

        //unpackaging (maybe needs to be reworked by using references in macroparam func)

        std::vector<T>& n_v   = params.n;
        std::vector<T>& u_1_v = params.u_1;
        std::vector<T>& T_v   = params.T;
        std::vector<T>& q_1_v = params.q_1;

        for (int x_i = 0; x_i < N_x; ++x_i) {
            T sqrt_T = std::sqrt(T_v[x_i]);
                // T nu = n_v[x_i] * sqrt_T * data.delta;  // for ball molecule model only!
                T nu = n_v[x_i] * T_v[x_i] * data.delta / mu(T_v[x_i]);

            T n_mult_inv_sqrt_Pi_T = n_v[x_i] * inv_sqrt_Pi / sqrt_T;
            T n_mult_sqrt_T_inv_Pi = n_v[x_i] * inv_sqrt_Pi * sqrt_T;

            for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
                T delta_xi_j = (xi_v[xi_j] - u_1_v[x_i]);
                T M_exp = std::exp(-delta_xi_j * delta_xi_j / T_v[x_i]);

                J_h[x_i][xi_j] = nu * (n_mult_inv_sqrt_Pi_T * M_exp - state.h[x_i][xi_j]);
                J_g[x_i][xi_j] = nu * (n_mult_sqrt_T_inv_Pi * M_exp - state.g[x_i][xi_j]);
            }
        }
        
        for (int xi_j = 0; xi_j < N_xi; ++xi_j) {
            xisRightSide(state.h, J_h, xi_j, diff_state.h, h_l[xi_j], h_r[xi_j]);
            xisRightSide(state.g, J_g, xi_j, diff_state.g, g_l[xi_j], g_r[xi_j]);
        }

        return diff_state;
    });

    T t{ 0.f };
    for (T t{ 0 }; data.t_end - t > 0; t += data.dt) {
        // Macroparameters1d<T> params = bgk1dMacroparameters(state, dxi, xi_v, apprInt);

        if (outputRule(t)) {
            output.print<T>(t, x_v, xi_v, bgk1dMacroparameters(state, dxi, xi_v, apprInt), state);
        }

        state = RkMethod::stepY(fullRightSide, state, std::min(data.dt, data.t_end - t));
    }

    output.print<T>(t, x_v, xi_v, bgk1dMacroparameters(state, dxi, xi_v, apprInt), state);
}
