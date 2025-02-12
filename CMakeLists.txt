cmake_minimum_required(VERSION 3.10)
project(HdbscanLibrary)

add_library(hdbscan STATIC
    src/HdbscanAlgorithm.cpp
    src/UndirectedGraph.cpp
    src/HdbscanCluster.cpp
)

target_include_directories(hdbscan PUBLIC include)

set_target_properties(hdbscan PROPERTIES CXX_STANDARD 17)

add_executable(test_hdbscan test/FourProminentClusterExample.cpp)
target_link_libraries(test_hdbscan PRIVATE hdbscan)