cmake_minimum_required(VERSION 3.22)
project(mathcpp LANGUAGES CXX)

add_library(mathcpp INTERFACE include/mathcpp/matrix.hpp
            include/mathcpp/vector.hpp)
target_include_directories(mathcpp INTERFACE include)
