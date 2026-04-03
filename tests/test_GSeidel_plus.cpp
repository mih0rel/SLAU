#include "gtest/gtest.h"
#include "solvers/GSeidel_plus.hpp"
#include "solvers/Chebyshev.hpp"


TEST(SYMMETRIC_GSEIDEL, TEST_1) {
    std::vector<std::vector<double>> a{{10, 1, 0}, {1,7, 0.1},{0,0,1}};
    CSR mat(a);
    std::vector<double> x{0,0,0};
    std::vector<double> b{20, 30, 1};
    std::vector<double>x_1 = symmetric_gseidel(mat,b, x, 1e-15);
    std::vector<double> ans{1.595652173913, 4.043478260869, 1.00000000000};
    std::cout << ans << std::endl;
    for (std::size_t i = 0; i < ans.size(); i ++) {
        ASSERT_NEAR(ans[i], x_1[i], 1e-8);
    }
    std::cout << "end of sim GS" << std::endl;
}

TEST(ACCELERATED_GSEIDEL, TEST_2) {
    std::vector<std::vector<double>> a{{10, 1, 0}, {1,7, 0.1},{0,0,1}};
    CSR mat(a);
    std::vector<double> x{0,0,0};
    std::vector<double> b{20, 30, 1};
    double lambda_max = calc_lambda_max(mat, 1e-8);
    std::cout << lambda_max << std::endl;
    std::vector<double>x_1 = accelerated_gseidel(mat,b, x, lambda_max * 2, 1e-15);
    std::cout << "i'm there" << std::endl;
    std::vector<double> ans{1.595652173913, 4.043478260869, 1.00000000000};
    for (std::size_t i = 0; i < ans.size(); i ++) {
        ASSERT_NEAR(ans[i], x_1[i], 1e-8);
    }
}
