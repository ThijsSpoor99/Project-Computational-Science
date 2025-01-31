# CMake configuration for the VTK example: solar_system_animateVTK
#
# This script sets up the build process for an executable that animates the 
# solar system using the VTK library. 
#
# Requirements:
# - VTK library must be installed.
# - Ensure the VTK lib and include directories are properly set in the system path.
#
# This CMake configuration:
# - Sets up a C++17 project.
# - Applies compiler-specific optimizations.
# - Finds and links the required VTK libraries.
#
# The workflow is as follows:
# - create a build directory and navigate to it.
# - run cmake .. to configure the project.
# - run cmake --build . to build the project. (some warnings may appear, but they can be ignored)
# - navigate to Debug/Release directory and run the executable.


cmake_minimum_required(VERSION 3.10)

# Set the project name and language
project(ProjectComputationalScience VERSION 1.0 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler-specific optimization flags
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /arch:AVX2 /Ob0 /fp:fast")
endif()

# Find and link the VTK package
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Add the executable
add_executable(solar_system_animateVTK solar_system_animateVTK.cpp)

# Link the VTK libraries
target_link_libraries(solar_system_animateVTK ${VTK_LIBRARIES})