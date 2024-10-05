cmake_minimum_required(VERSION 3.10)
project(Wavefunction)

# Manually set the C++17 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# Set Eigen include directory correctly to the 'eigen3' folder
set(libs_PATH "/groups/ngolubev/Packages/libs/")
set(EIGEN3_INCLUDE_DIR "${libs_PATH}/Eigen/eigen-3.4.0")
include_directories(${EIGEN3_INCLUDE_DIR})

# Set C++ compilation flags for performance optimizations
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -mavx512f")

# Enable OpenMP
find_package(OpenMP REQUIRED)

# Set source files
set(SOURCES_DENSITY density.cpp)
set(SOURCES_DENSITY_REDUCED density_reduced.cpp)

# Add the executable for extract.cpp
add_executable(density ${SOURCES_DENSITY})

# Add the executable for extract_pare.cpp
add_executable(density_reduced ${SOURCES_DENSITY_REDUCED})

# Link OpenMP to all executables and link stdc++fs for filesystem support
target_link_libraries(density OpenMP::OpenMP_CXX stdc++fs)
target_link_libraries(density_reduced OpenMP::OpenMP_CXX stdc++fs)

# Set runtime output directory for all executables
set_target_properties(density PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_target_properties(density_reduced PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Print confirmation messages
message(STATUS "Eigen include directory is set to: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "C++17 standard has been set using -std=c++17 flag.")
