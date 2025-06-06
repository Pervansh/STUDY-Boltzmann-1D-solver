#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <assert.h>
#include <concepts>
#include <cmath>

#include <omp.h>

#include "BoltzmannEq1DSolvers.h"
#include "Mkt.h"
#include "MacroparametersFileInterface.h"
#include "rkMethods.h"
#include "BoltzmannEq1DProblems.h"

void ompSetup() {
    omp_set_num_threads(20);
    std::cout << "num_threads: " << omp_get_max_threads() << std::endl;
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

void emittingWallTest() {
    double t_end = 10;
    double dt = 0.001;

    const auto testRule = [dt](double t) { return (int)std::round(t / dt) % 100 == 0; };

    Bgk1dProblemData<double> data = emittingWallProblem<double>(
        1., 1.,
        1.5, 2.,
        400, 12.5,
        50, 6.,
        dt, t_end,
        100
    );

    Full1dStateOutput outputBall("output\\bgk1d\\testEmittingWallBall");
    bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, outputBall, testRule, ballMoleculeViscosityRule);

    Full1dStateOutput outputMaxwell("output\\bgk1d\\testEmittingWallMaxwell");
    bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, outputMaxwell, testRule, maxwellMoleculeViscosityRule);
}

void evaporatingWallTest() {
    double t_end = 5.;
    double dt = 0.0001;
    //double t_end = 100000 * dt;

    int t_step_cnt = (int)std::round((t_end + 2 * dt) / dt);
    int control_t_step = t_step_cnt / 5;

    const auto testRule = [dt, control_t_step](double t) { return (int)std::round(t / dt) % control_t_step == 0; };

    Bgk1dProblemData<double> data = evaporatingWallProblem<double>(
        1., 1.,
        2., 10.,
        150, 15.,
        100, 15.,
        dt, t_end,
        knToDelta(0.01)
    );

    /*Full1dStateOutput outputBall("output\\bgk1d\\testEvaporatingWallBall");
    bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, outputBall, testRule, ballMoleculeViscosityRule);*/

    //Full1dStateOutput outputMaxwell("output\\bgk1d\\testEvaporatingWallMaxwell");
    // bgk1dMethod<ExplicitEulerRK>(data, simpsonInt, outputMaxwell, testRule, maxwellMoleculeViscosityRule);

    Bgk1dProblemData<double> data2 = evaporatingWallProblem<double>(
        1., 1.,
        2., 10.,
        500, 15.,
        100, 9.,
        dt, t_end,
        knToDelta(0.01)
    );

    //Full1dStateOutput outputMaxwell2("output\\bgk1d\\testEvaporatingWallMaxwell2");
    //bgk1dMethod<ExplicitEulerRK>(data2, simpsonInt, outputMaxwell2, testRule, maxwellMoleculeViscosityRule);

    //Full1dStateOutput outputMaxwell3("output\\bgk1d\\testEvaporatingWallMaxwell3");
    //bgk1dMethod<ExplicitEulerRK>(data2, Cell1stOrderInt, outputMaxwell3, testRule, maxwellMoleculeViscosityRule);

    Bgk1dProblemData<double> data4 = evaporatingWallProblem<double>(
        1., 1.,
        2., 10.,
        300, 15.,
        100, 9.,
        dt, t_end,
        knToDelta(0.01)
    );

    Full1dStateOutput outputMaxwell4("output\\bgk1d\\testEvaporatingWallMaxwell4");
    bgk1dMethod<ExplicitEulerRK>(data4, simpsonInt, outputMaxwell4, testRule, maxwellMoleculeViscosityRule);
}

void density10Test() {
    double t_end = 0.2 * std::sqrt(2.);
    double dt = 0.00002;

    int t_step_cnt = (int)std::round(t_end / dt);
    int control_t_step = t_step_cnt / 100;

    Bgk1dProblemData<double> data = densityRiemannProblem<double>(1., 1., 0.125, 150 * 5, 0., 1., 45 * 8, 8.0, dt, t_end, 10000.);

    Full1dStateOutput output("output\\bgk1d\\testDensity10Multiple");
    const auto testRule = [dt, control_t_step](double t) { return (int)std::round(t / dt) % control_t_step == 0; };
    // bgk1dMethod<ExplicitEulerRK>(data, Cell1stOrderInt, output, testRule);
    bgk1dMethod<ExplicitEulerRK>(data, simpsonInt, output, testRule);   
}

void sodTest() {
    double t_end = 0.2 * std::sqrt(2.);
    double dt = 0.00002;

    int t_step_cnt = (int)std::round(t_end / dt);
    int control_t_step = t_step_cnt / 100;

    // Bgk1dProblemData<double> data = densityTemperatureRiemannProblem<double>(1., 1., 0.8, 0.125, 150 * 5, 0., 1., 45 * 8, 8.0, dt, t_end, 10000.);
    Bgk1dProblemData<double> data = densityTemperatureRiemannProblem<double>(1., 1., 0.8, 0.125, 150 * 3, 0., 1., 45 * 4, 8.0, dt, t_end, 10000.);

    Full1dStateOutput output("output\\bgk1d\\SOD_1");
    const auto testRule = [dt, control_t_step](double t) { return (int)std::round(t / dt) % control_t_step == 0; };
    bgk1dMethod<ExplicitEulerRK>(data, simpsonInt, output, testRule);
}

int main() {
    ompSetup();

    //emittingWallTest();

    //density10Test();

    //evaporatingWallTest();

    sodTest();

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