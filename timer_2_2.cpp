#include "src/matrices/CSR.hpp"

#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cmath>

template<typename T>
double residual_norm(CSR<T>& matrix, const std::vector<double>& x, const std::vector<T>& b){
    std::vector<double> res = (matrix * x) - b;
    return std::sqrt(dot(res, res));
}


template<typename T>
std::vector<double> gseidel(CSR<T>& matrix, const std::vector<T>& b, const std::vector<T>& x_s, long double break_condition, int max_iter){
    std::ofstream file;
    file.open("GSeidel_2_2.csv");
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

    file.close();
    std::cout << "GSeidel iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}


template<typename T>
std::vector<double> jacobi(CSR<T>& matrix, const std::vector<T>& b, const std::vector<T>& x_s, long double break_condition, int max_iter){
    std::ofstream file;
    file.open("Jacobi_2_2.csv");
    file << "iter,residual,time_us" << std::endl;

    std::vector<double> x(x_s.size());
    double r_cur = 1e8;
    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        for(std::size_t j = 0; j < x_0.size(); j++){
            double tmp = 0;
            double d = 0;
            for(std::size_t p = matrix.get_rows(j); p < matrix.get_rows(j + 1); p++){
                if(matrix.get_cols(p) != j)
                    tmp += matrix.get_values(p) * x_0[matrix.get_cols(p)];
                else
                    d = matrix.get_values(p);
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
std::vector<double> symmetric_gseidel(CSR<T>& matrix, const std::vector<T>& b, const std::vector<T>& x_s, long double break_condition, int max_iter){
    std::ofstream file;
    file.open("SymmetricGSeidel_2_2.csv");
    file << "iter,residual,time_us" << std::endl;

    double r_cur = 1e8;
    int k = 0;
    std::vector<double> x_0(x_s.size());
    std::vector<double> x_old(x_s.size());
    std::copy(x_s.begin(), x_s.end(), x_0.begin());

    auto start = std::chrono::steady_clock::now();

    for(; r_cur > break_condition && k < max_iter;){
        std::copy(x_0.begin(), x_0.end(), x_old.begin());
        for(std::size_t i = 0; i < x_0.size(); i++){
            std::size_t j = x_0.size() - 1 - i;
            double L = 0;
            double U = 0;
            double d = 0;
            for(std::size_t p = matrix.get_rows(j); p < matrix.get_rows(j + 1); p++){
                if(matrix.get_cols(p) < j)
                    L += matrix.get_values(p) * x_old[matrix.get_cols(p)];
                else if(matrix.get_cols(p) > j)
                    U += matrix.get_values(p) * x_0[matrix.get_cols(p)];
                else
                    d = matrix.get_values(p);
            }
            x_0[j] = (b[j] - L - U) / d;
        }

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

    file.close();
    std::cout << "SymmetricGSeidel iterations = " << k << ", residual = " << r_cur << std::endl;
    return x_0;
}


int main(){
    std::vector<std::vector<double>> a{
        {18, 5, 2, 1, 0, 0, 0, 0},
        {5, 18, 5, 2, 1, 0, 0, 0},
        {2, 5, 18, 5, 2, 1, 0, 0},
        {1, 2, 5, 18, 5, 2, 1, 0},
        {0, 1, 2, 5, 18, 5, 2, 1},
        {0, 0, 1, 2, 5, 18, 5, 2},
        {0, 0, 0, 1, 2, 5, 18, 5},
        {0, 0, 0, 0, 1, 2, 5, 18}
    };
    CSR mat(a);
    std::vector<double> x{0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> b{1, 2, 3, 4, 5, 6, 7, 8};

    long double break_condition = 1e-15;
    int max_iter = 100000;

    std::vector<double> x_1 = jacobi(mat, b, x, break_condition, max_iter);
    std::vector<double> x_2 = gseidel(mat, b, x, break_condition, max_iter);
    std::vector<double> x_3 = symmetric_gseidel(mat, b, x, break_condition, max_iter);

    std::cout << x_1;
    std::cout << x_2;
    std::cout << x_3;

    return 0;
}
