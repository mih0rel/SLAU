#include "gtest/gtest.h"
#include "solvers/GSeidel.hpp"


TEST(gseidel, TEST_1) {
    std::vector<std::vector<double>> a{{12, 5,6}, {1,10, 8},{1,2,4}};
    CSR mat(a);
    std::vector<double> x{0,0,0};
    std::vector<double> b{7,4,1};
    std::vector<double>x_1 = gseidel(mat,b, x, 1e-15);
    std::vector<double> ans{0.4461538461538462, 0.4076923076923077, -0.0653846153846154};
    for (std::size_t i = 0; i < ans.size(); i ++) {
        ASSERT_NEAR(ans[i], x_1[i], 1e-8);
    }
}
