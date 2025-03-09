#pragma once

#include <vector>
#include <algorithm>
#include <assert.h>

template<class T>
std::vector<T> sum(const std::vector<T>& a, const std::vector<T>& b) {
    int n = std::min(a.size(), b.size());
    std::vector<T> res(n);

    for (int i = 0; i < n; ++i) {
        res[i] = a[i] + b[i];
    }

    return res;
}

template<class T>
std::vector<T> diff(const std::vector<T>& a, const std::vector<T>& b) {
    int n = std::min(a.size(), b.size());
    std::vector<T> res(n);

    for (int i = 0; i < n; ++i) {
        res[i] = a[i] - b[i];
    }

    return res;
}

template<class T, class U = T>
std::vector<T> div(const std::vector<T>& a, U coef) {
    int n = a.size();
    std::vector<T> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] / coef;
    }
    return res;
}

template<class T, class U = T>
std::vector<T> mul(const std::vector<T>&a, U coef) {
    int n = a.size();
    std::vector<T> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] * coef;
    }
    return res;
}

template<class T>
inline std::vector<T> operator+ (const std::vector<T>& a, const std::vector<T>& b) {
    return sum(a, b);
}

template<class T>
inline std::vector<T> operator- (const std::vector<T>& a, const std::vector<T>& b) {
    return diff(a, b);
}

template<class T>
std::vector<T> operator- (const std::vector<T>& a) {
    int n = a.size();
    std::vector<T> res(n);

    for (int i = 0; i < n; ++i) {
        res[i] = -a[i];
    }

    return res;
}

template<class T, class U = T>
inline std::vector<T> operator* (const std::vector<T>& v, U coef) {
    return mul(v, coef);
}

template<class T, class U = T>
inline std::vector<T> operator* (U coef, const std::vector<T>& v) {
    return mul(v, coef);
}

template<class T, class U = T>
inline std::vector<T> operator/ (const std::vector<T>& v, U coef) {
    return div(v, coef);
}

template<class T>
inline std::vector<T> operator+= (std::vector<T>& a, const std::vector<T>& b) {
    return a = sum(a, b);
}

template<class T>
inline std::vector<T> operator-= (std::vector<T>& a, const std::vector<T>& b) {
    return a = diff(a, b);
}

template<class T, class U = T>
inline std::vector<T> operator*= (std::vector<T>& v, U coef) {
    return v = mul(v, coef);
}

template<class T, class U = T>
inline std::vector<T> operator/= (std::vector<T>& v, U coef) {
    return v = div(v, coef);
}

template<typename T, typename Stream>
std::vector<T> readVector(Stream& stream, int size) {
    std::vector<T> readed(size);
    for (int i = 0; i < size; ++i) {
        stream >> readed[i];
    }

    return readed;
}

template<typename T, typename Stream>
void printVector(const std::vector<T>& v, Stream& stream, std::string sep = " ") {
    if (v.size() == 0) return;
    
    stream << v[0];
    for (int i = 0; i < v.size(); ++i) {
        stream << sep << v[i];
    }
}

template<typename T>
std::vector<T> rangeVector(int n, T a, T b) {
    assert(n > 1);

    std::vector<T> res(n);
    res[0] = a;

    T diff = (a - b) / n;

    for (int i = 1; i < n; ++i) {
        res[i] = res[i - 1] + diff;
    }

    return res;
}

template<typename T>
std::vector<T> centerRangeVector(int n, T a, T b) {
    assert(n > 0);

    std::vector<T> res(n);

    T diff = (a - b) / n;
    res[0] = a + diff / 2;

    for (int i = 1; i < n; ++i) {
        res[i] = res[i - 1] + diff;
    }

    return res;
}
