cmake_minimum_required(VERSION 4.3)
project(LeastSquaresApproximation)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

set(CMAKE_CXX_STANDARD 20)

add_executable(LeastSquaresApproximation main.cpp
        LeastSquaresApproximation.cpp
        LeastSquaresApproximation.h
        demos.cpp
        demos.h)
