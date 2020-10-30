cmake_minimum_required(VERSION 3.14)
project(K-mer)

set(MACOSX_RPATH ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

include_directories(src/include)
set(SOURCE_FILES
        src/BipartiteGraph.cpp
        src/Bipartition.cpp
        src/Geodesic.cpp
        src/PhyloTree.cpp
        src/PhyloTreeEdge.cpp
        src/Distance.cpp
        src/Ratio.cpp
        src/RatioSequence.cpp
        src/Tools.cpp)

add_executable(tests ${SOURCE_FILES} src/test.cpp src/bitset_hash.h)
add_executable(Kmers ${SOURCE_FILES} src/TestMain.cpp  src/bitset_hash.h headers/species.h headers/dist_matrix.h headers/neighbour_joining.h headers/phylogenetic_tree.h headers/utilities.h src/species.cpp src/dist_matrix.cpp src/neighbour_joining.cpp src/phylogenetic_tree.cpp src/utilities.cpp)
add_executable(build_tree ${SOURCE_FILES} src/build_tree.cpp src/bitset_hash.h)

#find_program(GENGETOPT_EXECUTABLE gengetopt)
#add_executable(Kmers main.cpp headers/species.h headers/dist_matrix.h headers/neighbour_joining.h headers/phylogenetic_tree.h headers/utilities.h src/species.cpp src/dist_matrix.cpp src/neighbour_joining.cpp src/phylogenetic_tree.cpp src/utilities.cpp)
#add_executable(Kmers  )
