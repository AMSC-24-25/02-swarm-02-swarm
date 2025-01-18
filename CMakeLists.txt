cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(02-Swarm VERSION 0.1)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(cmake/modules/CPM.cmake)

CPMAddPackage(
    NAME benchmark
    GITHUB_REPOSITORY google/benchmark
    VERSION 1.9.1
    OPTIONS "BENCHMARK_ENABLE_TESTING Off"
)

if(benchmark_ADDED)
    # enable c++11 to avoid compilation errors with Google Benchmark
    set_target_properties(benchmark PROPERTIES CXX_STANDARD 11)
endif()

CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    VERSION 1.15.2
    OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt"
)

file(GLOB_RECURSE SOURCES "src/*.cpp")

include_directories(include)

add_executable(main exec/main.cpp ${SOURCES})


enable_testing()

add_executable(prova exec/test.cpp)
target_link_libraries(prova PRIVATE gtest)

add_test(NAME RunTests COMMAND Tests)

set(BENCHMARK_SOURCES
    bench/benchmark_test.cpp
)
add_executable(bench ${BENCHMARK_SOURCES})
target_link_libraries(bench PRIVATE benchmark)

message(STATUS "")
message(STATUS "Project ${PROJECT_NAME} v${PROJECT_VERSION} successfully configured.")
message(STATUS "")
