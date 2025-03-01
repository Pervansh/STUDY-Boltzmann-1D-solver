#include <iostream>
#include <fstream>
#include <memory>

#include "AdvectionEq1DSolvers.h"
#include "HyperbolicFileInterface.h"
#include "MonotoneFluxes.h"
#include "Fluxes.h"

#define DEBUG_PRINT_TYPE

template <typename T, typename RkMethod, typename MonotoneFluxType, typename AdvectionFluxType>
void selectInterpolation(
    std::ostream& output,
    const Advection1dTaskFileData<T>& data,
    const MonotoneFluxType& monotoneFlux,
    const TaskData<T, AdvectionFluxType>& taskData
) {
    if (data.interpolationMethodId == 0) {
        uniformMinmodRecVecRkMethod<T, RkMethod>(taskData, monotoneFlux, output);
    } else if (data.interpolationMethodId == 1) {
        uniformConstRecVecRkMethod<T, RkMethod>(taskData, monotoneFlux, output);
    } else if (data.interpolationMethodId == 2) {
        std::clog << "[DEBUG]: WENO5 interpolation";
        uniformWeno5RecVecRkMethod<T, RkMethod>(taskData, monotoneFlux, output);
    } else {
        std::cerr << "[ERROR]: bad interpolationMethodId: " << data.interpolationMethodId << std::endl;
    }
}

template <typename T, typename RkMethod, typename MonotoneFluxType>
void selectFlux(
    std::ostream& output,
    const Advection1dTaskFileData<T>& data,
    const MonotoneFluxType& monotoneFlux)
{
#ifdef DEBUG_PRINT_TYPE
    std::clog << "[DEBUG]: RkMethod type: " << typeid(RkMethod).name() << '\n';
#endif // DEBUG_PRINT_TYPE

    if (data.fluxId == 0) {
        LinearFlux<T> flux(T(1.f));
        auto taskData = constructTaskData(data, flux);
        selectInterpolation<T, RkMethod>(output, data, monotoneFlux, taskData);

    } else {
        std::cerr << "[ERROR]: bad fluxId: " << data.fluxId << std::endl;
    }
}

template <typename T, typename RkMethod>
void selectMonotoneFlux(std::ostream& output, const Advection1dTaskFileData<T>& data){
    if (data.monotomeFluxId == 0) {
        // musclSelectRkMethod<T, decltype(rusanovFlux)>(data);
        auto ptr = rusanovFlux<T>;
        selectFlux<T, RkMethod>(output, data, ptr);

    } else if (data.monotomeFluxId == 1) {
        // musclSelectRkMethod<T, UniformSimpleLaxFriedrichsFlux>(data);
        UniformSimpleLaxFriedrichsFlux<T> laxFriedrichsFlux(data.dx, data.dt);
        selectFlux<T, RkMethod>(output, data, laxFriedrichsFlux);

    } else if (data.monotomeFluxId == 2) {
        // musclSelectRkMethod<T, UniformSimpleLaxFriedrichsFlux>(data);
        AccurateFlux<T> accurateFlux(true);
        selectFlux<T, RkMethod>(output, data, accurateFlux);

    } else {
        std::cerr << "[ERROR]: bad monotomeFluxId: " << data.monotomeFluxId << std::endl;
    }
}

template <typename T>
void selectRkMethod(std::ostream& output, const Advection1dTaskFileData<T>& data) {
    if (data.rkMethodId == 0) {
        selectMonotoneFlux<T, HeunsMethodRK>(output, data);

    } else if (data.rkMethodId == 1) {
        selectMonotoneFlux<T, SSPRK3>(output, data);

    } else if (data.rkMethodId == 2) {
        selectMonotoneFlux<T, ExplicitEulerRK>(output, data);

    } else {
        std::cerr << "[ERROR]: bad rkMethodId: " << data.rkMethodId << std::endl;
    }
}

template <typename T>
void launchScheme(const Advection1dTaskFileData<T>& data) {
    std::ofstream output(data.outputFileName + ".txt");

    selectRkMethod(output, data);
    
    output.close();
}

int main() {
    std::ifstream input("input.txt");
    auto data = advection1dTaskRead<double>(input);
    input.close();

    launchScheme<double>(data);

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