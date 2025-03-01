#pragma once

#include <vector>
#include <iostream>

template <typename T>
using rightSideFunc = vector<T>(*)(T, const vector<T>&);

template <typename T>
using RkVectorStepTY = vector<T>(*)(rightSideFunc<T>, T, const vector<T>&, T);

template<typename T>
inline std::vector<T> pendulumFunc(T, const std::vector<T>& xs) {
    auto x1 = xs[0];
    auto x2 = xs[1];

    return { x2 , -x1 };
}

template <typename T>
void testRkMethodTY(
    RkVectorStepTY<T> step, 
    rightSideFunc<T> func,
    std::istream& input,
    std::ostream& output
) {
    int n;
    input >> n;
    std::vector<T> y0s;
    for (auto& yi : y0s) input >> yi;
    int i 
}
