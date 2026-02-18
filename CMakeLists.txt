cmake_minimum_required(VERSION 3.20)
project(SLE)

enable_testing()

# Sources live in repository root; keep CMake simple and build here.
add_library(sle_lib
    TridiagonalMatrix.cpp
)
target_include_directories(sle_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(GTest QUIET)
if (GTest_FOUND)
    add_executable(test_TridiagonalMatrixSolver
        test_TridiagonalMatrixSolver.cpp
    )
    target_link_libraries(test_TridiagonalMatrixSolver PRIVATE sle_lib GTest::gtest_main)
    target_include_directories(test_TridiagonalMatrixSolver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

    add_test(NAME test_TridiagonalMatrixSolver COMMAND test_TridiagonalMatrixSolver)
endif()
