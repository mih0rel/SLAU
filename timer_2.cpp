#include "src/matrices/CSR.hpp"

#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cmath>
//Если timer.cpp был так хорош, где timer_2.cpp?
//(⊙_⊙)
template<typename T>
double residual_norm(CSR<T>& matrix, const std::vector<double>& x, const std::vector<T>& b){
    std::vector<double> res = (matrix * x) - b;
    return std::sqrt(dot(res, res));
}


template<typename T>
std::vector<double> gseidel(CSR<T>& matrix, const std::vector<T>& b, const std::vector<T>& x_s, long double break_condition, int max_iter){
    std::ofstream file;
    file.open("GSeidel.csv");
    file << "iter,residual,time_us" << std::endl;

    double tmp = 0;
    double r_cur = 1e8;
    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        for(std::size_t j = 0; j < x_0.size(); j++){
            tmp = 0;
            double d = 0;
            for(std::size_t p = matrix.get_rows(j); p < matrix.get_rows(j + 1); p++){
                if(matrix.get_cols(p) != j){
                    tmp += matrix.get_values(p) * x_0[matrix.get_cols(p)];
                }
                else{
                    d = matrix.get_values(p);
                }
            }
            x_0[j] = (b[j] - tmp) / d;
        }

        r_cur = residual_norm(matrix, x_0, b);
        auto finish = std::chrono::steady_clock::now();
        long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << r_cur << "," << time_us << std::endl;
    }

    file.close();
    std::cout << "GSeidel iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}


template<typename T>
std::vector<double> jacobi(CSR<T>& matrix, const std::vector<T>& b, const std::vector<T>& x_s, long double break_condition, int max_iter){
    std::ofstream file;
    file.open("Jacobi.csv");
    file << "iter,residual,time_us" << std::endl;

    std::vector<double> x(x_s.size());
    double tmp = 0;
    double r_cur = 1e8;
    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        for(std::size_t j = 0; j < x_0.size(); j++){
            tmp = 0;
            double d = 0;
            for(std::size_t p = matrix.get_rows(j); p < matrix.get_rows(j + 1); p++){
                if(matrix.get_cols(p) != j){
                    tmp += matrix.get_values(p) * x_0[matrix.get_cols(p)];
                }
                else{
                    d = matrix.get_values(p);
                }
            }
            x[j] = (b[j] - tmp) / d;
        }

        std::copy(x.begin(), x.end(), x_0.begin());
        r_cur = residual_norm(matrix, x_0, b);
        auto finish = std::chrono::steady_clock::now();
        long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << r_cur << "," << time_us << std::endl;
    }

    file.close();
    std::cout << "Jacobi iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}


template<typename T>
std::vector<double> FPI(CSR<T>& matrix, const std::vector<T>& b, const std::vector<T>& x_s, double tau, long double break_condition, int max_iter){
    std::ofstream file;
    file.open("SimpleIterations.csv");
    file << "iter,residual,time_us" << std::endl;

    std::vector<double> x(x_s.size());
    double r_cur = 1e8;
    int k = 0;
    std::vector<double> res(x_s.size());
    std::vector<double> x_0(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        res = (matrix * x_0) - b;
        x = x_0 - tau * res;
        std::copy(x.begin(), x.end(), x_0.begin());

        r_cur = std::sqrt(dot(res, res));
        auto finish = std::chrono::steady_clock::now();
        long long time_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

        k++;
        file << k << "," << r_cur << "," << time_us << std::endl;
    }

    file.close();
    std::cout << "SimpleIterations iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}


int main(){
    std::vector<std::vector<double>> a{
        {10, 2, 3, 1, 1},
        {9, 20, 1, 6, 1},
        {9, 7, 30, 2, 3},
        {3, 4, 10, 40, 12},
        {1, 3, 6, 2, 50}
    };
    CSR mat(a);
    std::vector<double> x{1, 1, 1, 1, 1};
    std::vector<double> b{4, 2, 2, 3, 8};

    long double break_condition = 1e-15;
    int max_iter = 100000;

    std::vector<double> x_1 = gseidel(mat, b, x, break_condition, max_iter);
    std::vector<double> x_2 = FPI(mat, b, x, 0.018, break_condition, max_iter);
    std::vector<double> x_3 = jacobi(mat, b, x, break_condition, max_iter);

    std::cout<< x_1;
    std::cout<< x_2;
    std::cout<<x_3;

    return 0;
}
