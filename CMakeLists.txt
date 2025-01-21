cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(02-Swarm VERSION 0.1)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_TESTS "Configures the test suite." ON)
option(ENABLE_BENCHMARKS "Configures the benchmark suite." ON)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Default build type." FORCE)
endif()

# Default flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "Default flags for Debug builds." FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "Default flags for Release builds." FORCE)

include(cmake/modules/CPM.cmake)

find_package(OpenMP REQUIRED)

file(GLOB_RECURSE SOURCES "src/*.cpp")

include_directories(include)

set(WARNINGS -Wall -Wextra -Wpedantic -Werror -Wshadow)

add_executable(main exec/main.cpp ${SOURCES})
target_compile_options(main PUBLIC ${WARNINGS})
target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)

if(ENABLE_TESTS)
    set(GOOGLE_TEST_VERSION "1.15.2" CACHE STRING "" FORCE)
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION ${GOOGLE_TEST_VERSION}
        OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt"
    )
    enable_testing()

    add_executable(test_function_convergence test/test_function_convergence.cpp ${SOURCES})
    target_compile_options(test_function_convergence PUBLIC ${WARNINGS})
    target_link_libraries(test_function_convergence PUBLIC OpenMP::OpenMP_CXX gtest gtest_main)
    include(GoogleTest)
    gtest_discover_tests(test_function_convergence)

    add_executable(TestRastrigin test/TestRastrigin.cpp ${SOURCES})
    target_compile_options(TestRastrigin PUBLIC ${WARNINGS})
    target_link_libraries(TestRastrigin PUBLIC OpenMP::OpenMP_CXX gtest gtest_main)

    add_executable(TestRosenbrock test/TestRosenbrock.cpp ${SOURCES})
    target_compile_options(TestRosenbrock PUBLIC ${WARNINGS})
    target_link_libraries(TestRosenbrock PUBLIC OpenMP::OpenMP_CXX gtest gtest_main)
endif()

if(ENABLE_BENCHMARKS)
    set(GOOGLE_BENCHMARK_VERSION "1.9.1" CACHE STRING "" FORCE)
    CPMAddPackage(
        NAME benchmark
        GITHUB_REPOSITORY google/benchmark
        VERSION ${GOOGLE_BENCHMARK_VERSION}
        OPTIONS "BENCHMARK_ENABLE_TESTING Off"
    )

    if(benchmark_ADDED)
        # enable c++11 to avoid compilation errors with Google Benchmark
        set_target_properties(benchmark PROPERTIES CXX_STANDARD 11)
    endif()

    set(BENCHMARK_SOURCES
        bench/bench_swarm_search.cpp
    )
    add_executable(bench ${BENCHMARK_SOURCES} ${SOURCES})
    target_link_libraries(bench PRIVATE benchmark)
endif()

# Clang Format configuration
find_program(CLANG_FORMAT_EXECUTABLE
    NAMES clang-format
    clang-format-14
    clang-format-15
    clang-format-16
    clang-format-17
    clang-format-18
    clang-format-19
    DOC "Path to the clang-format executable."
)

if(CLANG_FORMAT_EXECUTABLE)
    message(STATUS "Found Clang Format at ${CLANG_FORMAT_EXECUTABLE}. Adding new target 'format'.")
    file(GLOB_RECURSE ALL_SOURCE_FILES
        bench/*.cpp exec/*.cpp include/*.hpp src/*.cpp test/*.cpp)
    add_custom_target(format COMMAND ${CLANG_FORMAT_EXECUTABLE} --style=file -i ${ALL_SOURCE_FILES})
else()
    message(WARNING "Clang Format not found in PATH.")
endif()

message(STATUS "")
message(STATUS "Project ${PROJECT_NAME} v${PROJECT_VERSION} successfully configured.")
message(STATUS "CMAKE_BUILD_TYPE         : ${CMAKE_BUILD_TYPE}")
message(STATUS "ENABLE_TESTS             : ${ENABLE_TESTS}")

if(ENABLE_TESTS)
    message(STATUS "Google Test version      : ${GOOGLE_TEST_VERSION}")
endif()

message(STATUS "ENABLE_BENCHMARKS        : ${ENABLE_BENCHMARKS}")

if(ENABLE_BENCHMARKS)
    message(STATUS "Google Benchmark version : ${GOOGLE_BENCHMARK_VERSION}")
endif()

message(STATUS "")
