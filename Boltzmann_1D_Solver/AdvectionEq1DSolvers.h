#pragma once

#include <algorithm>
#include <cmath>
#include <ostream>
#include <memory>
#include <vector>

#include "AdditionalMath.h"
#include "rkMethods.h"
#include "Interpolation.h"
#include "VectorOperations.h"

template <typename T, typename AdvectionFluxType>
/*
    концепт для флакс тайпа: flux(T) -> T
*/
struct TaskData {
    // Left boundry of computation domain
    T a;
    // Right boundry of computation domain
    T b;
    // Spatial step
    T dx;
    // Time step
    T dt;
    // Time of simulation end
    T tEnd;
    // Flux functor
    const AdvectionFluxType& flux;
    // Number of elements in u0 array (i.e. a number of cells). Required: N >= 2!
    int N;
    // Initial data array (u0[i] = u0_{i}, where i - cell number), i = 0, ..., N - 1
    const T* u0;

    TaskData(T a, T b, T dx, T dt, T tEnd, const AdvectionFluxType& flux, int N, const T* u0)
        : a(a), b(b), dx(dx), dt(dt), tEnd(tEnd), flux(flux), N(N), u0(u0)
    {}
};

template <typename T, typename RkMethod, typename MonotoneFluxType, typename AdvectionFluxType>
void uniformMinmodRecRkMethod(
    const TaskData<T, AdvectionFluxType>& taskData,
    const MonotoneFluxType& monotoneFlux,
    std::ostream& output
) {
    const auto& N = taskData.N;

    std::unique_ptr<T[]> us = std::make_unique<T[]>(taskData.N);
    std::unique_ptr<T[]> usPrev = std::make_unique<T[]>(taskData.N);

    // cell centre printing
    for (int i = 0; i < N; ++i) {
        output << taskData.a + (2. * i + 1.) / 2. * taskData.dx << ' ';
    }
    output << '\n';

    // 0th layer printing
    output << 0. << ' ';
    for (int i = 0; i < N; ++i) {
        us[i] = taskData.u0[i];
        output << us[i] << ' ';
    }
    output << '\n';

    // Number of cells is N => number of boundries is (N + 1)
    // array of interpolated (extrapolated) values of u_i on right (i - 1)-th cell's border, uLs[i] = u^L_{i - 1/2}
    std::unique_ptr<T[]> uLs = std::make_unique<T[]>(taskData.N + 1);
    // array of interpolated (extrapolated) values of u_i on left i-th cell's border, uRs[i] = u^R_{i - 1/2}
    std::unique_ptr<T[]> uRs = std::make_unique<T[]>(taskData.N + 1);
    // array of values of monotone fluxes on cell's borders, fs[i] = f_{i - 1/2}
    std::unique_ptr<T[]> fs  = std::make_unique<T[]>(taskData.N + 1);

    T t = 0;
    for (; t <= taskData.tEnd; t += taskData.dt) {
        // Time layers changing
        std::swap(us, usPrev);

        /*
        // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
        T halfDelta0 = T(0.5f) * minmod(usPrev[0], usPrev[1] - usPrev[0]);
        T uL0 = u[i] + halfDelta; // Left value for right cell boundry
        T uR0 = u[i] - halfDelta; // Right value for left cell boundry

        fs[0] = monotoneFlux(taskData.flux(uL0), taskData.flux(uR0), uL0, uR0);


        // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
        T halfDelta0 = T(0.5f) * minmod(usPrev[N - 1] - usPrev[N - 2], -usPrev[N - 1]);
        T uLN = u[i] + halfDelta; // Left value for right cell boundry
        T uRN = u[i] - halfDelta; // Right value for left cell boundry

        fs[N - 1] = monotoneFlux(taskData.flux(uLN), taskData.flux(uRN), uLN, uRN);

        for (int i = 1; i < N - 1; ++i) {
            // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
            T halfDelta = T(0.5f) * minmod(usPrev[i] - usPrev[i - 1], usPrev[i + 1] - usPrev[i]);
            T uL = u[i] + halfDelta; // Left value for right cell boundry
            T uR = u[i] - halfDelta; // Right value for left cell boundry
            
            fs[i] = monotoneFlux(taskData.flux(uL), taskData.flux(uR), uL, uR);
        }
        */

        // Domain boundry conditions and treatment (periodic)

        // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
        T halfDelta0 = T(0.5f) * minmod(usPrev[0] - usPrev[N - 1], usPrev[1] - usPrev[0]); // periodic condition
        uLs[1] = usPrev[0] + halfDelta0; // Left value for right cell boundry
        uRs[0] = usPrev[0] - halfDelta0; // Right value for left cell boundry
        uRs[N] = uRs[0]; // periodic condition

        // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
        T halfDeltaN = T(0.5f) * minmod(usPrev[N - 1] - usPrev[N - 2], usPrev[0] - usPrev[N - 1]); // periodic condition
        uLs[N]     = usPrev[N - 1] + halfDeltaN; // Left value for right cell boundry
        uRs[N - 1] = usPrev[N - 1] - halfDeltaN; // Right value for left cell boundry
        uLs[0] = uLs[N]; // periodic condition

        // Calculations of uLs and uRs for inner boundries (iterating over cells)
        for (int i = 1; i < N - 1; ++i) {
            // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
            T halfDelta = T(0.5f) * minmod(usPrev[i] - usPrev[i - 1], usPrev[i + 1] - usPrev[i]);
            uLs[i + 1] = usPrev[i] + halfDelta; // Left value for right cell boundry
            uRs[i]     = usPrev[i] - halfDelta; // Right value for left cell boundry
        }

        // Monotone fluxes calculation for all cell's boundries
        for (int i = 0; i <= N; ++i) {
            fs[i] = monotoneFlux(taskData.flux(uLs[i]), taskData.flux(uRs[i]), uLs[i], uRs[i]);
        }

        int i = 0; // cell index for RK
        auto rkRightPart = [&i, &taskData, &fs](T u) {
            return (fs[i] - fs[i + 1]) / taskData.dx;
        };

        // TVD RK step to obtain next time layer solution us
        for (; i < N; ++i) {
            // us[i] = RkMethod<T, T>::stepY(rkRightPart, usPrev[i], taskData.dt);
            us[i] = RkMethod::stepY(rkRightPart, usPrev[i], taskData.dt);
        }
    }

    // last layer printing
    output << t - taskData.dt << ' ';
    for (int i = 0; i < N; i++) {
        output << us[i] << ' ';
    }
    output << '\n';
}

template <typename T, typename RkMethod, typename MonotoneFluxType, typename AdvectionFluxType>
void uniformConstRecRkMethod(
    const TaskData<T, AdvectionFluxType>& taskData,
    const MonotoneFluxType& monotoneFlux,
    std::ostream& output
) {
#ifdef DEBUG_PRINT_TYPE
    std::clog << "[DEBUG]: RkMethod type: " << typeid(RkMethod).name() << '\n';
#endif // DEBUG_PRINT_TYPE

    const auto& N = taskData.N;

    std::unique_ptr<T[]> us     = std::make_unique<T[]>(taskData.N);
    std::unique_ptr<T[]> usPrev = std::make_unique<T[]>(taskData.N);

    // cell centre printing
    for (int i = 0; i < N; ++i) {
        output << taskData.a + (2. * i + 1.) / 2. * taskData.dx << ' ';
    }
    output << '\n';

    // 0th layer printing
    output << 0. << ' ';
    for (int i = 0; i < N; ++i) {
        us[i] = taskData.u0[i];
        output << us[i] << ' ';
    }
    output << '\n';

    // Number of cells is N => number of boundries is (N + 1)
    // array of values of monotone fluxes on cell's borders, fs[i] = f_{i - 1/2}
    std::unique_ptr<T[]> fs = std::make_unique<T[]>(taskData.N + 1);

    T t = 0;
    for (; t <= taskData.tEnd; t += taskData.dt) {
        // Time layers changing
        std::swap(us, usPrev);
        
        // Domain boundry conditions treatment (periodic)
        fs[0]     = monotoneFlux(taskData.flux(usPrev[N - 1]), taskData.flux(usPrev[0]), usPrev[N - 1], usPrev[0]);
        fs[N] = monotoneFlux(taskData.flux(usPrev[0]), taskData.flux(usPrev[N - 1]), usPrev[0], usPrev[N - 1]);

        // Monotone fluxes calculation for all cell's boundries
        for (int i = 1; i < N; ++i) {
            fs[i] = monotoneFlux(taskData.flux(usPrev[i - 1]), taskData.flux(usPrev[i]), usPrev[i - 1], usPrev[i]);
        }

        int i = 0; // cell index for RK
        auto rkRightPart = [&i, &taskData, &fs](T u) {
            return (fs[i] - fs[i + 1]) / taskData.dx;
            };

        // TVD RK step to obtain next time layer solution us
        for (; i < N; ++i) {
            // us[i] = RkMethod<T, T>::stepY(rkRightPart, usPrev[i], taskData.dt);
            us[i] = RkMethod::stepY(rkRightPart, usPrev[i], taskData.dt);
        }
    }

    // last layer printing
    output << t - taskData.dt << ' ';
    for (int i = 0; i < N; i++) {
        output << us[i] << ' ';
    }
    output << '\n';
}

template <typename T, typename RkMethod, typename MonotoneFluxType, typename AdvectionFluxType>
void uniformConstRecVecRkMethod(
    const TaskData<T, AdvectionFluxType>& taskData,
    const MonotoneFluxType& monotoneFlux,
    std::ostream& output
) {
#ifdef DEBUG_PRINT_TYPE
    std::clog << "[DEBUG]: RkMethod type: " << typeid(RkMethod).name() << '\n';
#endif // DEBUG_PRINT_TYPE

    const auto& N = taskData.N;

    std::vector<T> us(N);
    std::vector<T> usPrev(N);

    // cell centre printing
    for (int i = 0; i < N; ++i) {
        output << taskData.a + (2. * i + 1.) / 2. * taskData.dx << ' ';
    }
    output << '\n';

    // 0th layer printing
    output << 0. << ' ';
    for (int i = 0; i < N; ++i) {
        us[i] = taskData.u0[i];
        output << us[i] << ' ';
    }
    output << '\n';

    T t = 0;
    for (; t <= taskData.tEnd; t += taskData.dt) {
        // Time layers changing
        std::swap(us, usPrev);

        // Operator, generating all right parts for all u_{i}
        auto rkRightPartOp = [&taskData, &N, &monotoneFlux](const std::vector<T>& us) {
            // array of values of monotone fluxes on cell's borders, fs[i] = f_{i - 1/2}
            std::vector<T> fs(N + 1); // Can be gotten rid of it, if code is optimised
            std::vector<T> diffs(N + 1);

            // Domain boundry conditions treatment (periodic)
            fs[0] = monotoneFlux(taskData.flux(us[N - 1]), taskData.flux(us[0]), us[N - 1], us[0]);
            fs[N] = monotoneFlux(taskData.flux(us[0]), taskData.flux(us[N - 1]), us[0], us[N - 1]);

            // Monotone fluxes calculation for all cell's boundries
            for (int i = 1; i < N; ++i) {
                fs[i] = monotoneFlux(taskData.flux(us[i - 1]), taskData.flux(us[i]), us[i - 1], us[i]);
            }

            for (int i = 0; i < N; ++i) {
                diffs[i] = (fs[i] - fs[i + 1]) / taskData.dx;
            }

            return std::move(diffs);
        };

        us = RkMethod::stepY(rkRightPartOp, usPrev, taskData.dt);
    }

    // last layer printing
    output << t - taskData.dt << ' ';
    for (int i = 0; i < N; i++) {
        output << us[i] << ' ';
    }
    output << '\n';
}

template <typename T, typename RkMethod, typename MonotoneFluxType, typename AdvectionFluxType>
void uniformMinmodRecVecRkMethod(
    const TaskData<T, AdvectionFluxType>& taskData,
    const MonotoneFluxType& monotoneFlux,
    std::ostream& output
) {
#ifdef DEBUG_PRINT_TYPE
    std::clog << "[DEBUG]: RkMethod type: " << typeid(RkMethod).name() << '\n';
#endif // DEBUG_PRINT_TYPE

    const auto& N = taskData.N;

    std::vector<T> us(N);
    std::vector<T> usPrev(N);

    // cell centre printing
    for (int i = 0; i < N; ++i) {
        output << taskData.a + (2. * i + 1.) / 2. * taskData.dx << ' ';
    }
    output << '\n';

    // 0th layer printing
    output << 0. << ' ';
    for (int i = 0; i < N; ++i) {
        us[i] = taskData.u0[i];
        output << us[i] << ' ';
    }
    output << '\n';

    // Number of cells is N => number of boundries is (N + 1)
    // array of interpolated (extrapolated) values of u_i on right (i - 1)-th cell's border, uLs[i] = u^L_{i - 1/2}
    std::vector<T> uLs(N + 1);
    // array of interpolated (extrapolated) values of u_i on left i-th cell's border, uRs[i] = u^R_{i - 1/2}
    std::vector<T> uRs(N + 1);
    // array of values of monotone fluxes on cell's borders, fs[i] = f_{i - 1/2}
    std::vector<T> fs(N + 1);

    T t = 0;
    for (; t <= taskData.tEnd; t += taskData.dt) {
        // Time layers changing
        std::swap(us, usPrev);

        // Operator, generating all right parts for all u_{i}
        auto rkRightPartOp = [&taskData, &uLs, &uRs, &fs, &N, &monotoneFlux](const std::vector<T>& usPrev) {
            std::vector<T> diffs(N + 1);

            // Domain boundry conditions and treatment (periodic)

            // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
            T halfDelta0 = T(0.5f) * minmod(usPrev[0] - usPrev[N - 1], usPrev[1] - usPrev[0]); // periodic condition
            uLs[1] = usPrev[0] + halfDelta0; // Left value for right cell boundry
            uRs[0] = usPrev[0] - halfDelta0; // Right value for left cell boundry
            uRs[N] = uRs[0]; // periodic condition

            // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
            T halfDeltaN = T(0.5f) * minmod(usPrev[N - 1] - usPrev[N - 2], usPrev[0] - usPrev[N - 1]); // periodic condition
            uLs[N] = usPrev[N - 1] + halfDeltaN; // Left value for right cell boundry
            uRs[N - 1] = usPrev[N - 1] - halfDeltaN; // Right value for left cell boundry
            uLs[0] = uLs[N]; // periodic condition

            // Calculations of uLs and uRs for inner boundries (iterating over cells)
            for (int i = 1; i < N - 1; ++i) {
                // half delta is without dx^{-1} because final formulae for uL & uR do not have dx
                T halfDelta = T(0.5f) * minmod(usPrev[i] - usPrev[i - 1], usPrev[i + 1] - usPrev[i]);
                uLs[i + 1] = usPrev[i] + halfDelta; // Left value for right cell boundry
                uRs[i] = usPrev[i] - halfDelta; // Right value for left cell boundry
            }

            // Monotone fluxes calculation for all cell's boundries
            for (int i = 0; i <= N; ++i) {
                fs[i] = monotoneFlux(taskData.flux(uLs[i]), taskData.flux(uRs[i]), uLs[i], uRs[i]);
            }

            for (int i = 0; i < N; ++i) {
                diffs[i] = (fs[i] - fs[i + 1]) / taskData.dx;
            }

            return std::move(diffs);
        };

        us = RkMethod::stepY(rkRightPartOp, usPrev, taskData.dt);
    }

    // last layer printing
    output << t - taskData.dt << ' ';
    for (int i = 0; i < N; i++) {
        output << us[i] << ' ';
    }
    output << '\n';
}

template <typename T, typename RkMethod, typename MonotoneFluxType, typename AdvectionFluxType>
void uniformWeno5RecVecRkMethod(
    const TaskData<T, AdvectionFluxType>& taskData,
    const MonotoneFluxType& monotoneFlux,
    std::ostream& output
) {
#ifdef DEBUG_PRINT_TYPE
    std::clog << "[DEBUG]: RkMethod type: " << typeid(RkMethod).name() << '\n';
#endif // DEBUG_PRINT_TYPE

    const auto& N = taskData.N;

    /*
        Looped array declaration!
        If array ys[N] of y_i is looped, then it means
        ys[i + 2] = y_i, y[0] = y_{N - 2}, y[1] = y_{N - 1}, y[N + 2] = y_0, y[N + 3] = y_1
    */

    std::vector<T> us(N);
    std::vector<T> usPrev(N);
    std::vector<T> loopedUsPrev(N + 4);

    // Looped array of x_i
    std::vector<T> xs(N + 4);

    // cell centre printing
    for (int i = 0; i < N; ++i) {
        xs[i + 2] = (2. * i + 1.) / 2. * taskData.dx;
        output << taskData.a + (2. * i + 1.) / 2. * taskData.dx << ' ';
    }
    output << '\n';

    // Looping xs
    xs[0]     = xs[N - 2];
    xs[1]     = xs[N - 1];
    xs[N + 2] = xs[0];
    xs[N + 3] = xs[1];

    // 0th layer printing
    output << 0. << ' ';
    for (int i = 0; i < N; ++i) {
        us[i] = taskData.u0[i];
        output << us[i] << ' ';
    }
    output << '\n';

    // Number of cells is N => number of boundries is (N + 1)
    // array of interpolated (extrapolated) values w/ WENO5 of u_i on right (i - 1)-th cell's border, uLs[i] = u^L_{i - 1/2}
    std::vector<T> uLs(N + 1);
    // array of interpolated (extrapolated) values w/ WENO5 of u_i on left i-th cell's border, uRs[i] = u^R_{i - 1/2}
    std::vector<T> uRs(N + 1);
    // array of values of monotone fluxes on cell's borders, fs[i] = f_{i - 1/2}
    std::vector<T> fs(N + 1);
    // array of interpolated (extrapolated) values of u_i on left i-th cell's border, vLs[i][r] = u^{(r), L}_{i - 1/2}
    std::vector<std::vector<T>> vLs(N + 1, std::vector<T>(3));
    // array of interpolated (extrapolated) values of u_i on left i-th cell's border, vRs[i][r] = u^{(r), R}_{i - 1/2}
    std::vector<std::vector<T>> vRs(N + 1, std::vector<T>(3));

    // Operator, generating all right parts for all u_{i}
    auto rkRightPartOp = [&taskData, &uLs, &uRs, &fs, &N, &xs, &loopedUsPrev, &monotoneFlux, &vLs, &vRs](
        const std::vector<T>& usPrev
    ) {
        std::vector<T> diffs(N + 1); // array for right parts for current RK step

        for (int i = 0; i < N; ++i) {
            loopedUsPrev[i + 2] = usPrev[i];
        }

        // Looping loopedUsPrev
        loopedUsPrev[0] = loopedUsPrev[N - 2];
        loopedUsPrev[1] = loopedUsPrev[N - 1];
        loopedUsPrev[N + 2] = loopedUsPrev[0];
        loopedUsPrev[N + 3] = loopedUsPrev[1];
        
        // Calculations of vLs and vRs for all boundries (iterating over cells)
        for (int i = 0; i < N; ++i) {
            // Note: not optimised, identical calculations in each interPolynom call
            for (int r = 0; r < 3; ++r) {
                // Left value for right cell boundry
                vLs[i + 1][r] = interPolynom3AtPoint(
                    xs.data() + i + 2 - r,
                    loopedUsPrev.data() + i + 2 - r,
                    T(i + 1) * taskData.dx);
                // Right value for left cell boundry
                vRs[i][r] = interPolynom3AtPoint(
                    xs.data() + i + 2 - r,
                    loopedUsPrev.data() + i + 2 - r,
                    T(i) * taskData.dx);
            }

            T ds[3] = { T(3.f) / 10, T(3.f) / 5, T(1.f) / 10 };
            T eps = 1e-6;
            T betas[3] = {
                T(13.f) / 12 * sqr(loopedUsPrev[i + 2] - 2 * loopedUsPrev[i + 3] + loopedUsPrev[i + 4]) +
                    T(0.25f) * sqr(3 * loopedUsPrev[i + 2] - 4 * loopedUsPrev[i + 3] + loopedUsPrev[i + 4]),
                T(13.f) / 12 * sqr(loopedUsPrev[i + 1] - 2 * loopedUsPrev[i + 2] + loopedUsPrev[i + 3]) +
                    T(0.25f) * sqr(loopedUsPrev[i + 1] + loopedUsPrev[i + 3]),
                T(13.f) / 12 * sqr(loopedUsPrev[i] - 2 * loopedUsPrev[i + 1] + loopedUsPrev[i + 2]) +
                    T(0.25f) * sqr(loopedUsPrev[i] - 4 * loopedUsPrev[i + 1] + 3 * loopedUsPrev[i + 2])
            };
            T alphas[3] = { ds[0] / sqr(eps + betas[0]), ds[1] / sqr(eps + betas[1]), ds[2] / sqr(eps + betas[2]) };
            T alphaSum = alphas[0] + alphas[1] + alphas[2];
            T omegas[3] = { alphas[0] / alphaSum, alphas[1] / alphaSum, alphas[2] / alphaSum };

            uLs[i + 1] = omegas[0] * vLs[i + 1][0] + omegas[1] * vLs[i + 1][1] + omegas[2] * vLs[i + 1][2];
            uRs[i] = omegas[0] * vRs[i][0] + omegas[1] * vRs[i][1] + omegas[2] * vRs[i][2];
        }

        // Domain boundry condition treatment (periodic)
        uLs[0] = uLs[N];
        uRs[N] = uRs[0];

        // Monotone fluxes calculation for all cell's boundries
        for (int i = 0; i <= N; ++i) {
             fs[i] = monotoneFlux(taskData.flux(uLs[i]), taskData.flux(uRs[i]), uLs[i], uRs[i]);
        }

        for (int i = 0; i < N; ++i) {
            diffs[i] = (fs[i] - fs[i + 1]) / taskData.dx;
        }

        return std::move(diffs);
    };

    T t = 0;
    ///*
    for (; t <= taskData.tEnd; t += taskData.dt) {
        // Time layers changing
        std::swap(us, usPrev);

        us = RkMethod::stepY(rkRightPartOp, usPrev, std::min(taskData.dt, taskData.tEnd - t));
    }
    //*/

    // last layer printing
    output << t - taskData.dt << ' ';
    for (int i = 0; i < N; i++) {
        output << us[i] << ' ';
    }
    output << '\n';
}
