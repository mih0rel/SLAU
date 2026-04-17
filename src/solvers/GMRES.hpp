#pragma once

#include "matrices/CSR.hpp"
#include "matrices/DenseMatrix.hpp"

#include <cmath>
#include <utility>
#include <vector>

template <typename T>
void ArnoldiIter(CSR<T> &matrix, DenseMatrix<T> &K,
                      DenseMatrix<T> &Hessenberg, const unsigned int i){
    std::vector<T> t = matrix * K.get_col(i);
    T h;
    for (std::size_t k = 0; k <= i; k++){
        h = dot(t, K.get_col(k));
        t = t - h * K.get_col(k);
        Hessenberg.change_elem(k, i) = h;
    }
    h = norm(t);
    Hessenberg.change_elem(i + 1, i) = h;
    t = 1 / h * t;
    K.write_col(t, i + 1);
}

template <typename T>
void GivensRotation(DenseMatrix<T> &Hessenberg,
                        std::vector<std::pair<double, double>>& rotations, const unsigned int i){
    for (std::size_t k = 0; k < i; k++){
        double h = rotations[k].first * Hessenberg(k, i) - rotations[k].second * Hessenberg(k + 1, i);
        double h_n = rotations[k].first * Hessenberg(k + 1, i) + rotations[k].second * Hessenberg(k, i);
        Hessenberg.change_elem(k, i) = h;
        Hessenberg.change_elem(k + 1, i) = h_n;
    }

    double cos_phi = Hessenberg(i, i) / std::sqrt(Hessenberg(i, i) * Hessenberg(i, i) + Hessenberg(i + 1, i) * Hessenberg(i + 1, i));
    double sin_phi = -Hessenberg(i + 1, i) / std::sqrt(Hessenberg(i, i) * Hessenberg(i, i) + Hessenberg(i + 1, i) * Hessenberg(i + 1, i));
    double h = cos_phi * Hessenberg(i, i) - sin_phi * Hessenberg(i + 1, i);
    double h_n = sin_phi * Hessenberg(i, i) + cos_phi * Hessenberg(i + 1, i);
    Hessenberg.change_elem(i, i) = h;
    Hessenberg.change_elem(i + 1, i) = h_n;

    rotations[i] = std::make_pair(cos_phi, sin_phi);
}

template<typename T>
std::vector<T> Gauss_Backward(const DenseMatrix<T>& matrix, std::vector<T> b) {
    std::vector<T> ans(b.size(), 0);
    ans[b.size() - 1] = b[b.size() - 1] / matrix(b.size() - 1, b.size() - 1);
    for (std::size_t i = b.size() - 1; i > 0; i--) {
        double tmp = b[i - 1];
        for (std::size_t k = b.size() - 1; k > i - 1; k--) tmp -= matrix(i - 1, k) * ans[k];
        ans[i - 1] = (1 / matrix(i - 1, i - 1)) * tmp;
    }
    return ans;
}

template<typename T>
std::pair<std::vector<double>, double> GMRES_Iter(CSR<T>& matrix, const std::vector<T>& b, std::vector<T> x_0, const unsigned int m){

    DenseMatrix<double> Hessenberg(0, m + 1, m);
    DenseMatrix<double> K(0, x_0.size(), m + 1);

    std::vector<std::pair<double, double>> rotations(m);

    std::vector<double> r_0 = matrix * x_0 - b;
    std::vector<double> e(m + 1, 0);
    e[0] = norm(r_0);

    K.write_col(1 / norm(r_0) * r_0, 0);

    for(unsigned int i = 0; i < m ; i++){
        ArnoldiIter(matrix, K ,Hessenberg, i);
        GivensRotation(Hessenberg, rotations, i);

        double e_p = rotations[i].first * e[i] - rotations[i].second * e[i + 1];
        double e_n = rotations[i].second * e[i] + rotations[i].first * e[i + 1];

        e[i] = e_p;
        e[i + 1] = e_n;
    }

    auto tmp = std::vector<T>(e.begin(), e.begin() + m);
    std::vector<double> y = Gauss_Backward(Hessenberg, tmp);
    std::vector<double> ans(x_0.size(), 0);
    for (unsigned int k = 0; k < y.size(); k++) ans = ans + y[k] * K.get_col(k);

    return {x_0 - ans, std::abs(e[m])};
}

template<typename T>
std::vector<double> GMRES(CSR<T>& matrix, const std::vector<T>& b,  const std::vector<T>& x_0, const unsigned int m, long double break_condition){
    std::pair<std::vector<double>, double> ans = GMRES_Iter(matrix, b, x_0, m);
    for(; ans.second > break_condition;) ans = GMRES_Iter(matrix, b, ans.first, m);
    return ans.first;
}
