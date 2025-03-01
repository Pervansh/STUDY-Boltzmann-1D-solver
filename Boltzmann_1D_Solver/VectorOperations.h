#pragma once

#include <vector>
#include <algorithm>

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

template<class T>
std::vector<T> div(const std::vector<T>& a, T coef) {
    int n = a.size();
    std::vector<T> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] / coef;
    }
    return res;
}

template<class T>
std::vector<T> mul(const std::vector<T>&a, T coef) {
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

template<class T>
inline std::vector<T> operator* (const std::vector<T>& v, T coef) {
    return mul(v, coef);
}

template<class T>
inline std::vector<T> operator* (T coef, const std::vector<T>& v) {
    return mul(v, coef);
}

template<class T>
inline std::vector<T> operator/ (const std::vector<T>& v, T coef) {
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

template<class T>
inline std::vector<T> operator*= (std::vector<T>& v, T coef) {
    return v = mul(v, coef);
}

template<class T>
inline std::vector<T> operator/= (std::vector<T>& v, T coef) {
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
