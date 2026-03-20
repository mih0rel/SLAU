#include "matrices/CSR.hpp"
#include <cmath>
#include <iostream>
#include <numbers>

template<typename T>
double calc_lambda_max(CSR<T>& matrix, double break_condition){
    std::cout << "i'm there CALC LM";
    double mu_cur = 0;
    double mu_prev = 1;
    std::vector<double> r(matrix.get_width(), 1);
    for(; std::abs(mu_cur - mu_prev) > break_condition;){
        mu_prev = mu_cur;
        r = 1 / norm(matrix * r) * matrix * r;
        mu_cur = dot(r, matrix * r) / dot(r, r);
    }
    return mu_cur;
}

std::vector<size_t> roots_order(int r){
    size_t n = static_cast<size_t>(std::pow(2, r));
    std::vector<size_t> indexes(n);
    indexes[static_cast<size_t>(std::pow(2, r - 1))] = 1;
    size_t step = 0;
    for (int i = 2; i <= r; ++i) {
        step = static_cast<size_t>(std::pow(2, r - i));
        for (size_t j = 0; j < n; j += 2 * step) {
            indexes[j + step] = static_cast<size_t>(std::pow(2, i)) - indexes[j] - 1;
        }
    }
    return indexes;
}

std::vector<double> calc_roots(int r, double lambda_min, double lambda_max){
    size_t n = static_cast<size_t>(std::pow(2, r));
    std::vector<double> roots(n);
    const double cos_n = std::cos(std::numbers::pi / static_cast<int>(n));
    const double sin_n = std::sin(std::numbers::pi/ static_cast<int>(n));
    const double cos_2n = std::cos(std::numbers::pi / (2 * static_cast<int>(n)));
    double sin_i = std::sin(std::numbers::pi/ (2 * static_cast<int>(n)));
    roots[0] = cos_2n;
    for(size_t i = 1; i < n / 2 + 1; i ++){
        roots[i] = roots[i - 1] * cos_n - sin_i * sin_n;
        sin_i = sin_i * cos_n + roots[i - 1] * sin_n;
        roots[n - i] = -roots[i - 1];
        roots[i - 1] = (lambda_min + lambda_max) / 2 + (lambda_max - lambda_min) / 2 * roots[i - 1];
        roots[n - i] = (lambda_min + lambda_max) / 2 + (lambda_max - lambda_min) / 2 * roots[n - i];
    }
    return roots;

}

template<typename T>
std::vector<double> Chebyshev_accel(CSR<T>& matrix, const std::vector<T>& b,  const std::vector<T>& x_s, int r, double lambda_min, double lambda_max, long double break_condition){
    std::vector<double> x(x_s.size());
    double r_cur = 1;
    std::vector<double> res(x_s.size()); 
    std::vector<double> x_0(x_s.size()) ;
    std::copy(x_s.begin(), x_s.end(), x_0.begin());
    std::vector<size_t> indexes = roots_order(r);
    std::vector<double> roots = calc_roots(r, lambda_min, lambda_max);
    for(; r_cur > break_condition;){
        for(auto& id:indexes) {
            res = (matrix * x_0 - b);
            x = x_0 - 1 / roots[id] * res;
            std::copy(x.begin(), x.end(), x_0.begin());
            r_cur = std::sqrt(dot(res, res)); 
            if (r_cur < break_condition) break;
        }
        r_cur = std::sqrt(dot(res, res));
    }
    return x_0;
}
