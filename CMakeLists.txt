project(mathcpp LANGUAGES CXX)

add_library(mathcpp INTERFACE include/mathcpp/matrix.hpp
            include/mathcpp/vector.hpp include/mathcpp/rk4.hpp)
target_include_directories(mathcpp INTERFACE include)
