#include "src/matrices/EllipticGenerator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <vector>

template<typename T>
double residual_norm(CSR<T>& matrix, const std::vector<double>& x, const std::vector<T>& b){
    std::vector<double> res = (matrix * x) - b;
    return std::sqrt(dot(res, res));
}

template<typename T>
double calc_lambda_max(CSR<T>& matrix, double break_condition){
    double mu_cur = 0;
    double mu_prev = 1;
    std::vector<double> r(matrix.get_width(), 1);

    for(; std::abs(mu_cur - mu_prev) > break_condition;){
        mu_prev = mu_cur;
        r = 1 / norm(matrix * r) * (matrix * r);
        mu_cur = dot(r, matrix * r) / dot(r, r);
    }

    return mu_cur;
}

double poisson_sor_parameter(std::size_t n){
    double mu = std::cos(std::numbers::pi / (n + 1));
    return 1 + std::pow(mu / (1 + std::sqrt(1 - mu * mu)), 2);
}

template<typename T>
std::vector<double> gseidel_iter(CSR<T>& matrix, const std::vector<T>& b, std::vector<double> x_0){
    std::vector<double> x_old(x_0.size());
    std::copy(x_0.begin(), x_0.end(), x_old.begin());

    for(std::size_t i = 0; i < x_0.size(); i++){
        std::size_t j = x_0.size() - 1 - i;
        double L = 0;
        double U = 0;
        double d = 0;
        for(std::size_t k = matrix.get_rows(j); k < matrix.get_rows(j + 1); k++){
            if(matrix.get_cols(k) < j)
                L += matrix.get_values(k) * x_old[matrix.get_cols(k)];
            else if(matrix.get_cols(k) > j)
                U += matrix.get_values(k) * x_0[matrix.get_cols(k)];
            else
                d = matrix.get_values(k);
        }
        x_0[j] = (b[j] - L - U) / d;
    }

    std::copy(x_0.begin(), x_0.end(), x_old.begin());
    for(std::size_t j = 0; j < x_0.size(); j++){
        double L = 0;
        double U = 0;
        double d = 0;
        for(std::size_t k = matrix.get_rows(j); k < matrix.get_rows(j + 1); k++){
            if(matrix.get_cols(k) < j)
                L += matrix.get_values(k) * x_0[matrix.get_cols(k)];
            else if(matrix.get_cols(k) > j)
                U += matrix.get_values(k) * x_old[matrix.get_cols(k)];
            else
                d = matrix.get_values(k);
        }
        x_0[j] = (b[j] - L - U) / d;
    }

    return x_0;
}

template<typename T>
std::vector<double> gseidel(CSR<T>& matrix, const std::vector<T>& b, const std::vector<double>& x_s, double break_condition, int max_iter){
    std::ofstream file("GSeidel_3.csv");
    file << "iter,residual,time_us" << std::endl;

    double r_cur = 1e8;
    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::vector<double> x_old(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        std::copy(x_0.begin(), x_0.end(), x_old.begin());
        for(std::size_t j = 0; j < x_0.size(); j++){
            double L = 0;
            double U = 0;
            double d = 0;
            for(std::size_t p = matrix.get_rows(j); p < matrix.get_rows(j + 1); p++){
                if(matrix.get_cols(p) < j)
                    L += matrix.get_values(p) * x_0[matrix.get_cols(p)];
                else if(matrix.get_cols(p) > j)
                    U += matrix.get_values(p) * x_old[matrix.get_cols(p)];
                else
                    d = matrix.get_values(p);
            }
            x_0[j] = (b[j] - L - U) / d;
        }

        r_cur = residual_norm(matrix, x_0, b);
        auto finish = std::chrono::steady_clock::now();
        long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << r_cur << "," << time_us << std::endl;
    }

    std::cout << "GSeidel iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}

template<typename T>
std::vector<double> sor(CSR<T>& matrix, const std::vector<T>& b, const std::vector<double>& x_s, double omega, double break_condition, int max_iter){
    std::ofstream file("SOR_3.csv");
    file << "iter,residual,time_us" << std::endl;

    double r_cur = 1e8;
    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::vector<double> x_old(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        std::copy(x_0.begin(), x_0.end(), x_old.begin());
        for(std::size_t j = 0; j < x_0.size(); j++){
            double L = 0;
            double U = 0;
            double d = 0;
            for(std::size_t p = matrix.get_rows(j); p < matrix.get_rows(j + 1); p++){
                if(matrix.get_cols(p) < j)
                    L += matrix.get_values(p) * x_0[matrix.get_cols(p)];
                else if(matrix.get_cols(p) > j)
                    U += matrix.get_values(p) * x_old[matrix.get_cols(p)];
                else
                    d = matrix.get_values(p);
            }
            x_0[j] = (1 - omega) * x_old[j] + omega * (b[j] - L - U) / d;
        }

        r_cur = residual_norm(matrix, x_0, b);
        auto finish = std::chrono::steady_clock::now();
        long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << r_cur << "," << time_us << std::endl;
    }

    std::cout << "SOR iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}

template<typename T>
std::vector<double> accelerated_gseidel(CSR<T>& matrix, const std::vector<T>& b, const std::vector<double>& x_s, double rho, double break_condition, int max_iter){
    std::ofstream file("ChebyshevGSeidel_3.csv");
    file << "iter,residual,time_us" << std::endl;

    std::vector<double> y_0(x_s.size());
    std::copy(x_s.begin(), x_s.end(), y_0.begin());

    double mu_0 = 1;
    double mu_1 = 1 / rho;
    int k = 1;

    auto start = std::chrono::steady_clock::now();

    std::vector<double> y_1 = gseidel_iter(matrix, b, y_0);
    double r_cur = residual_norm(matrix, y_1, b);
    auto finish = std::chrono::steady_clock::now();
    long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    file << k << "," << r_cur << "," << time_us << std::endl;

    std::vector<double> y(y_0.size());
    for(; r_cur > break_condition && k < max_iter;){
        y_1 = gseidel_iter(matrix, b, y_1);
        y = 2 * mu_1 / rho * y_1 - mu_0 * y_0;
        mu_0 = 2 / rho * mu_1 - mu_0;
        y /= mu_0;

        y_0 = y_1;
        y_1 = y;

        std::swap(mu_0, mu_1);

        r_cur = residual_norm(matrix, y_1, b);
        finish = std::chrono::steady_clock::now();
        time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << r_cur << "," << time_us << std::endl;
    }

    std::cout << "ChebyshevGSeidel iterations = " << k << ", residual = " << r_cur << std::endl;
    return y_1;
}

template<typename T>
std::vector<double> fastest_grad(CSR<T>& matrix, const std::vector<T>& b, const std::vector<double>& x_s, double break_condition, int max_iter){
    std::ofstream file("FastestGrad_3.csv");
    file << "iter,residual,time_us" << std::endl;

    double r_cur = 1e8;
    double tau = 0;
    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        std::vector<double> res = matrix * x_0 - b;
        tau = dot(res, res) / dot(res, matrix * res);
        x_0 = x_0 - tau * res;

        r_cur = residual_norm(matrix, x_0, b);
        auto finish = std::chrono::steady_clock::now();
        long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << r_cur << "," << time_us << std::endl;
    }

    std::cout << "FastestGrad iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}

template<typename T>
std::vector<double> cg(CSR<T>& matrix, const std::vector<T>& b, const std::vector<double>& x_s, double break_condition, int max_iter){
    std::ofstream file("CG_3.csv");
    file << "iter,residual,time_us" << std::endl;

    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    std::vector<double> r_cur = matrix * x_0 - b;
    std::vector<double> d = r_cur;

    auto start = std::chrono::steady_clock::now();

    for(; norm(r_cur) > break_condition && k < max_iter;){
        x_0 = x_0 - dot(r_cur, r_cur) / dot(d, matrix * d) * d;

        std::vector<double> r_prev = r_cur;
        r_cur = matrix * x_0 - b;
        d = r_cur + dot(r_cur, r_cur) / dot(r_prev, r_prev) * d;

        auto finish = std::chrono::steady_clock::now();
        long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << norm(r_cur) << "," << time_us << std::endl;
    }

    std::cout << "CG iterations = " << k << ", residual = " << norm(r_cur) << std::endl;
    return x_0;
}

int main(){
    std::size_t n = 20;
    double L = 1.0;
    auto f = [](double x, double y){return 1.0;};
    auto g = [](double x, double y){return 0.0;};

    Elliptic<double> problem = dirichlet<double>(n, L, f, g);
    CSR<double> mat = problem.matrix;
    std::vector<double> b = problem.rhs;
    std::vector<double> x(mat.get_width(), 0);

    double break_condition = 1e-8;
    int max_iter = 20000;
    double omega = poisson_sor_parameter(n);
    double lambda_max = calc_lambda_max(mat, 1e-8);
    double rho = 2 * lambda_max;

    std::cout << "omega = " << omega << std::endl;
    std::cout << "lambda_max = " << lambda_max << std::endl;
    std::cout << "rho = " << rho << std::endl;

    std::vector<double> x_1 = gseidel(mat, b, x, break_condition, max_iter);
    std::vector<double> x_2 = sor(mat, b, x, omega, break_condition, max_iter);
    std::vector<double> x_3 = accelerated_gseidel(mat, b, x, rho, break_condition, max_iter);
    std::vector<double> x_4 = fastest_grad(mat, b, x, break_condition, max_iter);
    std::vector<double> x_5 = cg(mat, b, x, break_condition, max_iter);

    std::cout << x_1;
    std::cout << x_2;
    std::cout << x_3;
    std::cout << x_4;
    std::cout << x_5;

    return 0;
}
