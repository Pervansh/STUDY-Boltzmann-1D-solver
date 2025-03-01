#pragma once

#include <iostream>
#include <memory>
#include <utility>
#include <string>

#include "AdvectionEq1DSolvers.h"

template <typename T>
struct Advection1dTaskFileData {
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
    // Number of elements in u0 array (i.e. a number of cells). required: N >= 2!
    int N;
    // Initial data array (u0[i] = u0_{i}, where i - cell number), i = 0, ..., N - 1
    std::shared_ptr<T[]> u0;

    int fluxId;
    int interpolationMethodId;
    int monotomeFluxId;
    int rkMethodId;

    std::string outputFileName;
};

template <typename T, typename AdvectionFluxType>
TaskData<T, AdvectionFluxType> constructTaskData(const Advection1dTaskFileData<T>& data, const AdvectionFluxType& flux) {
    return TaskData<T, AdvectionFluxType>(
        data.a, data.b, data.dx, data.dt, data.tEnd, flux, data.N, data.u0.get()
    );
}

template <typename T>
Advection1dTaskFileData<T> advection1dTaskRead(std::istream& input) {
    Advection1dTaskFileData<T> data;
    input >> data.a >> data.b >> data.dx >> data.dt >> data.tEnd >> data.N;

    data.u0 = std::make_shared<T[]>(data.N);
    for (int i = 0; i < data.N; ++i) {
        input >> data.u0[i];
    }

    input >> data.fluxId >> data.interpolationMethodId >> data.monotomeFluxId >> data.rkMethodId;
    input >> data.outputFileName;

    return std::move(data);
}
