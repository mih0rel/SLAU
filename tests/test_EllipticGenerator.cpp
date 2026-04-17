#include "gtest/gtest.h"
#include "matrices/EllipticGenerator.hpp"


TEST(ELLIPTIC_GENERATOR, TEST_1) {
    auto f = [](double x, double y){return 0.0;};
    auto g = [](double x, double y){return 0.0;};

    auto problem = dirichlet<double>(1, 1.0, f, g);

    ASSERT_EQ(problem.matrix.get_height(), 1);
    ASSERT_EQ(problem.matrix.get_width(), 1);
    ASSERT_NEAR(problem.h, 0.5, 1e-12);
    ASSERT_NEAR(problem.matrix(0, 0), 16.0, 1e-12);
    ASSERT_NEAR(problem.rhs[0], 0.0, 1e-12);
}

TEST(ELLIPTIC_GENERATOR, TEST_2) {
    auto f = [](double x, double y){return 1.0;};
    auto g = [](double x, double y){return 1.0;};

    auto problem = dirichlet<double>(2, 1.0, f, g);

    ASSERT_EQ(problem.matrix.get_height(), 4);
    ASSERT_EQ(problem.matrix.get_width(), 4);
    ASSERT_NEAR(problem.h, 1.0 / 3.0, 1e-12);

    ASSERT_NEAR(problem.matrix(0, 0), 36.0, 1e-12);
    ASSERT_NEAR(problem.matrix(0, 1), -9.0, 1e-12);
    ASSERT_NEAR(problem.matrix(0, 2), -9.0, 1e-12);
    ASSERT_NEAR(problem.matrix(0, 3), 0.0, 1e-12);

    for(std::size_t i = 0; i < problem.rhs.size(); i++){
        ASSERT_NEAR(problem.rhs[i], 19.0, 1e-12);
    }
}
