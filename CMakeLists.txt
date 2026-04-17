cmake_minimum_required(VERSION 3.20)
project(SLE LANGUAGES CXX)

enable_testing()

add_library(sle_lib
    src/matrices/TridiagonalMatrix.cpp
)
target_sources(sle_lib PUBLIC
    src/decompositions/Householder.hpp
    src/matrices/CSR.hpp
    src/matrices/DenseMatrix.hpp
    src/matrices/EllipticGenerator.hpp
    src/matrices/TridiagonalMatrix.hpp
    src/matrices/Vector.hpp
    src/solvers/GSeidel.hpp
    src/solvers/GSeidel_plus.hpp
    src/solvers/GaussQR.hpp
    src/solvers/GMRES.hpp
    src/solvers/Jacobi.hpp
    src/solvers/CG.hpp
    src/solvers/Chebyshev.hpp
    src/solvers/Fastest_Grad.hpp
    src/solvers/SOR.hpp
    src/solvers/SimpleIterations.hpp
    src/solvers/TridiagonalMatrixSolver.hpp
)
target_include_directories(sle_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_features(sle_lib PUBLIC cxx_std_20)
set_target_properties(sle_lib PROPERTIES CXX_EXTENSIONS OFF)

find_package(GTest QUIET)
if (GTest_FOUND)
    function(add_sle_test target_name source_name)
        add_executable(${target_name}
            ${source_name}
        )
        target_link_libraries(${target_name} PRIVATE sle_lib GTest::gtest_main)
        add_test(NAME ${target_name} COMMAND ${target_name})
    endfunction()

    add_sle_test(test_DenseMatrix tests/test_DenseMatrix.cpp)
    add_sle_test(test_CSR tests/test_CSR.cpp)
    add_sle_test(test_EllipticGenerator tests/test_EllipticGenerator.cpp)
    add_sle_test(test_TridiagonalMatrixSolver tests/test_TridiagonalMatrixSolver.cpp)
    add_sle_test(test_Jacobi tests/test_Jacobi.cpp)
    add_sle_test(test_Chebyshev tests/test_Chebyshev.cpp)
    add_sle_test(test_Fastest_Grad tests/test_Fastest_Grad.cpp)
    add_sle_test(test_GMRES tests/test_GMRES.cpp)
    add_sle_test(test_CG tests/test_CG.cpp)
    add_sle_test(test_SOR tests/test_SOR.cpp)
    add_sle_test(test_SimpleIterations tests/test_SimpleIterations.cpp)
    add_sle_test(test_GSeidel tests/test_GSeidel.cpp)
    add_sle_test(test_GSeidel_plus tests/test_GSeidel_plus.cpp)
else()
    message(STATUS "GTest not found, tests will not be built")
endif()
