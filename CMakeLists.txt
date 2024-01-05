cmake_minimum_required(VERSION 3.16)
project(genome_variant_simulator)

set(CMAKE_CXX_STANDARD 17)

add_executable(genome_variant_simulator main.cpp)

target_compile_options(genome_variant_simulator
        PUBLIC
        -Wall -Wextra -pedantic
        -O3 -funroll-loops -fomit-frame-pointer -ffast-math
)
