#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <assert.h>
#include <concepts>
#include <cmath>

#include "BoltzmannEq1DSolvers.h"
#include "Mkt.h"
#include "MacroparametersFileInterface.h"
#include "rkMethods.h"

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

void studentSpringDensityTest() {
    double t_end = 5;
    double dt = 0.001;
    Bgk1dProblemData<double> data = densityRiemannProblem<double>(1., 2., 1., 150, 5., 45, 4.0, dt, t_end, 100.);

    Full1dStateOutput output("output\\bgk1d\\testDensityMultiple");
    const auto testRule = [dt](double t) { return (int)std::round(t / dt) % 100 == 0; };
    bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, output, testRule);
}

void freeMoleculeDensityTest() {
    double t_end = 10;
    double dt = 0.001;

    const auto testRule = [dt](double t) { return (int)std::round(t / dt) % 100 == 0; };

    /*
    Bgk1dProblemData<double> data = densityRiemannProblem<double>(1., 2., 1., 150, 10., 46, 4.0, dt, t_end, 0.);
    Full1dStateOutput output("output\\bgk1d\\testFreeMoleculeDensity");
    bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, output, testRule);
    */

    ///*
    Bgk1dProblemData<double> shifted_data = densityRiemannProblem<double>(1., 2., 1., 150, 5., 224, 4.0, dt, t_end, 0., 1);
    Full1dStateOutput shifted_output("output\\bgk1d\\testFreeMoleculeDensityShifted");
    bgk1dMethod<ExplicitEulerRK>(shifted_data, Cell1stOrderInt, shifted_output, testRule);
    //*/
}

int main() {
    freeMoleculeDensityTest();

    /*
    Bgk1dProblemData<double> data = temperatureRiemannProblem<double>(1., 0.1, 6, 1., 5, 1.0, 0.01, 1., 1.);

    Macroparameters1dOutput output("output\\bgk1d\\testSmall");
    */

    //Bgk1dProblemData<double> data = temperatureRiemannProblem<double>(1., 0.1, 50, 1., 45, 4.0, 0.001, 0.1, 1.);

    //Macroparameters1dOutput output("output\\bgk1d\\testCrash");

    //Bgk1dProblemData<double> data = temperatureRiemannProblem<double>(1., 1., 50, 1., 45, 4., 0.01, 0.5, 1.);

    //Macroparameters1dOutput output("output\\bgk1d\\testConst");

    // t = 0.244
    // Bgk1dProblemData<double> data = densityRiemannProblem<double>(1., 2., 1., 50, 1., 45, 4.0, 0.001, 0.244, 20.);

    // Macroparameters1dOutput output("output\\bgk1d\\testDensity");
    // Full1dStateOutput output("output\\bgk1d\\testDensityMultiple");

    // const auto testRule = [](double t) { return t > 0.23; };

    // bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, output, testRule);

    //Full1dStateOutput output("output\\bgk1d\\testDensity");
    //bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, output);

    /*
    std::cout << "Hello world!" << std::endl;
    
    double a = 1.;
    auto linFlux { [&a](double u) {
        return a * u;
    } };

    std::unique_ptr<double[]> u0 = std::make_unique<double[]>(10);
    
    TaskData<double, decltype(linFlux)> taskData(0., 1., 0.1, 0.1, 1., linFlux, u0.get(), 10);

    float aF = 1.;
    auto linFluxF{ [&aF](float u) {
        return aF * u;
    } };

    std::unique_ptr<float[]> u0F = std::make_unique<float[]>(10);

    TaskData<float, decltype(linFluxF)> taskDataF(0.f, 1.f, 0.1f, 0.1f, 1.f, linFluxF, u0F.get(), 10);
    
    
    uniformMinmodRecRkMethod<float, SSPRK3>(
        taskDataF,
        rusanovFlux<float>,
        std::cout
    );
    */

    /*
    uniformMinmodRecRkMethod<double, HeunsMethodRK>(
        taskData,
        rusanovFlux<double>,
        std::cout
    );

    UniformSimpleLaxFriedrichsFlux<double> laxFriedrichsFlux(taskData.dx, taskData.dt);
    uniformMinmodRecRkMethod<double, SSPRK3>(
        taskData,
        laxFriedrichsFlux,
        std::cout
    );

    uniformMinmodRecRkMethod<double, HeunsMethodRK>(
        taskData,
        laxFriedrichsFlux,
        std::cout
    );
    */
    
    return 0;
}