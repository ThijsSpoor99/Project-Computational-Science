cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ProjectComputationalScience)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set optimization flags for Release build
# for release build do cmake -DCMAKE_BUILD_TYPE=Release ..
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
endif()

# Set debug symbols for Debug build
# for debug build do cmake -DCMAKE_BUILD_TYPE=Debug ..
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi")
endif()

# Find the VTK package
find_package(VTK REQUIRED)

find_package(OpenMP REQUIRED)

# Include the VTK headers
include(${VTK_USE_FILE})

# Add the executable
add_executable(solar_system_animateVTK solar_system_animateVTK.cpp)

# Link the VTK libraries
target_link_libraries(solar_system_animateVTK ${VTK_LIBRARIES})
target_link_libraries(solar_system_animateVTK OpenMP::OpenMP_CXX)