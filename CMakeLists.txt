cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(02-Swarm VERSION 0.1)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_TESTS "Configures the test suite." ON)
option(ENABLE_BENCHMARKS "Configures the benchmark suite." ON)
option(USE_MPI "Enables usage of MPI algorithms." OFF)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Default build type." FORCE)
endif()

# Default flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "Default flags for Debug builds." FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "Default flags for Release builds." FORCE)

include(cmake/modules/CPM.cmake)

find_package(OpenMP REQUIRED)

if(USE_MPI)
    find_package(MPI REQUIRED)
endif()

set(WARNINGS -Wall -Wextra -Wpedantic -Werror -Wshadow)

set(SOURCES
        src/Algorithm.cpp
        src/Creature.cpp
        src/EuclideanDistance.cpp
        src/GeneticAlgorithm.cpp
        src/Particle.cpp
        src/Rastrigin.cpp
        src/Rosenbrock.cpp
        src/Sphere.cpp
        src/Swarm.cpp
        src/MultiStochasticTunnelling.cpp
        src/Position.cpp
        src/DifferentialEvolution.cpp
        src/Candidate.cpp
        src/SimulatedAnnealing.cpp
        src/State.cpp
)
if(USE_MPI)
    list(APPEND SOURCES src/DistributedGeneticAlgorithm.cpp)
    list(APPEND SOURCES src/DistributedMultiStochasticTunnelling.cpp)
    list(APPEND SOURCES src/DistributedDifferentialEvolution.cpp)
    list(APPEND SOURCES src/DistributedSimulatedAnnealing.cpp)

endif()

# Define main library target
add_library(02-swarm ${SOURCES})
target_include_directories(02-swarm PUBLIC include)
target_compile_options(02-swarm PUBLIC ${WARNINGS})
target_link_libraries(02-swarm PUBLIC OpenMP::OpenMP_CXX)
if(USE_MPI)
    target_compile_definitions(02-swarm PUBLIC USE_MPI=1)
    target_link_libraries(02-swarm PUBLIC MPI::MPI_CXX)
endif()

add_executable(main exec/main.cpp)
if(USE_MPI)
    target_compile_definitions(main PUBLIC USE_MPI=1)
    target_link_libraries(main PUBLIC MPI::MPI_CXX)
endif()
target_link_libraries(main PUBLIC 02-swarm)

if(ENABLE_TESTS)
    set(GOOGLE_TEST_VERSION "1.15.2" CACHE STRING "" FORCE)
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION ${GOOGLE_TEST_VERSION}
        OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt"
    )
    enable_testing()

    include(GoogleTest)

    macro(add_test target)
        set(test_filename test_${target})
        add_executable(${test_filename} test/${test_filename}.cpp)
        target_link_libraries(${test_filename} PUBLIC 02-swarm gtest gtest_main)
        gtest_discover_tests(${test_filename})
    endmacro()

    add_test(swarm_convergence)
    add_test(genetic_convergence_omp)
    add_test(multi_tunnelling_convergence)
    add_test(de_convergence)
    add_test(sa_convergence)


    if(USE_MPI)
        add_executable(test_genetic_convergence_mpi test/test_genetic_convergence_mpi.cpp)
        target_compile_definitions(test_genetic_convergence_mpi PUBLIC USE_MPI=1)
        target_link_libraries(test_genetic_convergence_mpi PUBLIC 02-swarm gtest gtest_main MPI::MPI_CXX)
        gtest_discover_tests(test_genetic_convergence_mpi)

        add_executable(test_tunnelling_convergence_mpi test/test_multi_tunnelling_convergence_mpi.cpp test/main_test_de_mpi.cpp)
        target_compile_definitions(test_tunnelling_convergence_mpi PUBLIC USE_MPI=1)
        target_link_libraries(test_tunnelling_convergence_mpi PUBLIC 02-swarm gtest gtest_main MPI::MPI_CXX)
        gtest_discover_tests(test_tunnelling_convergence_mpi)

        add_executable(test_de_convergence_mpi
                test/test_de_convergence_mpi.cpp
                test/main_test_de_mpi.cpp
        )
        target_compile_definitions(test_de_convergence_mpi PUBLIC USE_MPI=1)
        target_link_libraries(test_de_convergence_mpi PUBLIC 02-swarm gtest MPI::MPI_CXX)
        gtest_discover_tests(test_de_convergence_mpi)

        add_executable(test_sa_convergence_mpi
                test/test_sa_convergence_mpi.cpp
                test/main_test_sa_mpi.cpp
        )
        target_compile_definitions(test_sa_convergence_mpi PUBLIC USE_MPI=1)
        target_link_libraries(test_sa_convergence_mpi PUBLIC 02-swarm gtest MPI::MPI_CXX)
        gtest_discover_tests(test_sa_convergence_mpi)

    endif(USE_MPI)

    add_executable(TestRastrigin test/TestRastrigin.cpp)
    target_link_libraries(TestRastrigin PUBLIC 02-swarm gtest gtest_main)

    add_executable(TestRosenbrock test/TestRosenbrock.cpp)
    target_link_libraries(TestRosenbrock PUBLIC 02-swarm gtest gtest_main)
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
        bench/bench_genetic_omp.cpp
        bench/bench_tunnelling_omp.cpp
        bench/bench_de_omp.cpp
        bench/bench_sa_omp.cpp

    )
    add_executable(bench ${BENCHMARK_SOURCES})
    target_link_libraries(bench PUBLIC 02-swarm benchmark benchmark_main)
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

message(STATUS "USE_MPI                  : ${USE_MPI}")

if(USE_MPI)
    message(STATUS "MPI version              : v${MPI_CXX_VERSION}")
    message(STATUS "MPIEXEC_EXECUTABLE       : ${MPIEXEC_EXECUTABLE}")
endif()

message(STATUS "")

