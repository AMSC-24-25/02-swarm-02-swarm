cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(Example VERSION 0.1)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(FetchContent)

# Disable tests
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_BENCHMARKS OFF CACHE BOOL "" FORCE)

# Download the library
FetchContent_Declare(
    02-Swarm
    GIT_REPOSITORY https://github.com/AMSC-24-25/02-swarm-02-swarm.git
    GIT_TAG main
)

# Configure the library
FetchContent_MakeAvailable(02-Swarm)

add_executable(example example.cpp)
target_link_libraries(example PUBLIC 02-swarm)
