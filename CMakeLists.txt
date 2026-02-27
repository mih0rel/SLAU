cmake_minimum_required(VERSION 3.20)
project(SLE)

enable_testing()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(sle_lib
    src/matrices/TridiagonalMatrix.cpp
)
target_include_directories(sle_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

find_package(GTest QUIET)
if (GTest_FOUND)
    add_executable(test_DenseMatrix
        tests/test_DenseMatrix.cpp
    )
    target_link_libraries(test_DenseMatrix PRIVATE GTest::gtest_main)
    target_include_directories(test_DenseMatrix PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_DenseMatrix COMMAND test_DenseMatrix)

    add_executable(test_CSR
        tests/test_CSR.cpp
    )
    target_link_libraries(test_CSR PRIVATE GTest::gtest_main)
    target_include_directories(test_CSR PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_CSR COMMAND test_CSR)

    add_executable(test_TridiagonalMatrixSolver
        tests/test_TridiagonalMatrixSolver.cpp
    )
    target_link_libraries(test_TridiagonalMatrixSolver PRIVATE sle_lib GTest::gtest_main)
    target_include_directories(test_TridiagonalMatrixSolver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME test_TridiagonalMatrixSolver COMMAND test_TridiagonalMatrixSolver)
endif()
