#include "gtest/gtest.h"
#include "solvers/CG.hpp"


TEST(CG, TEST_1) {
    std::vector<std::vector<double>> a{{10, 16, 4}, {17,150, 30},{3,30,300}};
    CSR mat(a);
    std::vector<double> x{1,1,1};
    std::vector<double> b{1,2,3};
 
    std::vector<double>x_1 = CG(mat, b, x, 1e-15);
    std::vector<double> ans{0.09521435692921244, 0.0007477567298105596, 0.008973080757726818};
    for (std::size_t i = 0; i < ans.size(); i ++) {
        ASSERT_NEAR(ans[i], x_1[i], 1e-8);
    }
}
