cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ProjectComputationalScience)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the VTK package
find_package(VTK REQUIRED)

# Include the VTK headers
include(${VTK_USE_FILE})

# Add the executable
add_executable(solar_system_animateVTK solar_system_animateVTK.cpp)

# Link the VTK libraries
target_link_libraries(solar_system_animateVTK ${VTK_LIBRARIES})