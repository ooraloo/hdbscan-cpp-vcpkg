cmake_minimum_required(VERSION 3.10)
project(HdbscanLibrary)

add_library(hdbscan STATIC
    src/HdbscanAlgorithm.cpp
    src/UndirectedGraph.cpp
    src/HdbscanCluster.cpp
)

target_include_directories(hdbscan PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(hdbscan PROPERTIES CXX_STANDARD 17)

add_executable(test_hdbscan test/FourProminentClusterExample.cpp)
target_link_libraries(test_hdbscan PRIVATE hdbscan)

install(TARGETS hdbscan
    EXPORT HDBSCANTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

install(EXPORT HDBSCANTargets
    FILE HDBSCANConfig.cmake
    NAMESPACE HDBSCAN::
    DESTINATION share/hdbscan
)