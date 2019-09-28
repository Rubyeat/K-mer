cmake_minimum_required(VERSION 3.14)
project(K-mer)

set(CMAKE_CXX_STANDARD 17)
find_program(GENGETOPT_EXECUTABLE gengetopt)
add_executable(Kmers main.cpp headers/species.h headers/dist_matrix.h headers/neighbour_joining.h headers/phylogenetic_tree.h headers/utilities.h src/species.cpp src/dist_matrix.cpp src/neighbour_joining.cpp src/phylogenetic_tree.cpp src/utilities.cpp src/rf.cpp)