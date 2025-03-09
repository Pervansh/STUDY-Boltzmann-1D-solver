#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <assert.h>
#include <concepts>

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

    data.a = x_max;
    data.b = -x_max;
    
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

int main() {
    Bgk1dProblemData<double> data = temperatureRiemannProblem<double>(1., 0.1, 50, 1., 50, 5., 0.01, 1., 1.);

    std::ofstream file("output\\test\\test_file.txt");
    file << "TEST!";
    file.close();

    Macroparameters1dOutput output("output\\bgk1d\\test");

    bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, output);

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