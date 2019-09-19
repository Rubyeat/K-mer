cmake_minimum_required(VERSION 3.14)
project(K-mer)

set(CMAKE_CXX_STANDARD 17)

add_executable(Kmers main.cpp headers/species.h src/species.cpp)