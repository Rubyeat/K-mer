cmake_minimum_required(VERSION 3.14)
project(K-mer)

set(CMAKE_CXX_STANDARD 17)

add_executable(Kmers main.cpp headers/species.h headers/dist_matrix.h headers/neighbour_joining.h headers/phylogenetic_tree.h headers/utilities.h src/species.cpp src/dist_matrix.c src/neighbour_joining.c src/phylogenetic_tree.c src/utilities.c)