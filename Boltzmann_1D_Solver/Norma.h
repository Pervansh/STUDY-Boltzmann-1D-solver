#pragma once

#include <vector>
#include <utility>
#include <cmath>
//#include "QuadMatrix.h"

// Т®ѓл дг≠™ж®© ≠Ѓађ §Ђп Ґ•™вЃаЃҐ ® ђ ва®ж бЃЃвҐ•вбвҐ•≠≠Ѓ

template<typename T>
using vectorNorm = T (*)(const std::vector<T>&);

template<class T>
class QuadMatrix;

template<typename T>
using matrixNorm = T (*)(const QuadMatrix<T>&);

// убическа  норма дл  вектора и матрицы
template <class T>
T norm_inf(const std::vector<T>& x) {
    size_t n = x.size();
    T xMax = 0;

    for (int i = 0; i < n; ++i) {
        if (fabs(x[i]) > xMax) {
            xMax = fabs(x[i]);
        }
    }
    return xMax;
}

/*
template <class T>
T norm_inf(const QuadMatrix<T>& A) {
    size_t n = A.order();
    T  aSum, aMax = 0;

    for (int i = 0; i < n; ++i) {
        aSum = 0;
        for (int j = 0; j < n; j++) {
            aSum += fabs(A(i, j));
        }
        if (aSum > aMax) {
            aMax = aSum;
        }
    }
    return aMax;
}
*/

//ќктаэдральна  норма дл  вектора и матрицы
template<class T>
T norm_1(const std::vector<T>& x) {
    size_t n = x.size();
    T xSum = 0;
    for (int i = 0; i < n; ++i) {
        xSum += fabs(x[i]);
    }
    return xSum;
}

/*
template<class T>
T norm_1(const QuadMatrix<T>& A) {
    size_t n = A.order();
    T aSum, aMax = 0;
    for (int j = 0; j < n; ++j) {
        aSum = 0;
        for (int i = 0; i < n; ++i) {
            aSum += fabs(A(i, j));
        }
        if (aSum > aMax) {
            aMax = aSum;
        }
    }
    return aMax;
}
*/

//Ўарова  норма дл  вектора и матрицы(норма ‘робениуса)
//template<class T>
//T norm_2(vector<T> x) {
//	size_t n = x.size();
//	T xSum = 0;
//	for (int i = 0; i < n; ++i) {
//		xSum += x[i] * x[i];
//	}
//	return sqrt(xSum);
//}
//Ёто не точно
//template<class T>
//T norm_F(QuadMatrix<T> A) {
//	size_t n = A.order();
//	T aSum = 0;
//	for (int i = 0; i < n; ++i) {
//		for (int j = 0; j < n; ++j) {
//			aSum += A(i, j) * A(i, j);
//		}
//	}
//	return sqrt(aSum);
//}
