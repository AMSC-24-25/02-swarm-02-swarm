# Project for AMSC 2024/2025 course
This is the project for the group 02. It implements some algorithms for which you can find a brief description in the following files:
 - [Swarm Search](https://github.com/AMSC-24-25/02-swarm-02-swarm/blob/main/SwarmSearch.md)
 - [Genetic Algorithm](https://github.com/AMSC-24-25/02-swarm-02-swarm/blob/main/GeneticAlgorithm.md)

For a detailed description of the Genetic Algorithm implementations, there is [a report](https://github.com/AMSC-24-25/02-swarm-02-swarm/blob/main/doc/report.pdf) available in the `doc` folder.

## How to compile and run
This project uses CMake as build system and requires at least version 3.14. You can install CMake by downloading it from [the official website](https://cmake.org/download).
You can compile it with:
```bash
cmake -S . -B build
cmake --build build -- all
```

Run with:
```bash
./build/main
```

Check that everything works correctly:
```bash
ctest --test-dir build
```

You can also run the benchmark suite by executing:
```bash
./build/bench
```

GoogleTest and Google Benchmark are automatically downloaded inside the build directory using [CPM](https://github.com/cpm-cmake/CPM.cmake), if they are enabled.

By default, the test suite and the benchmark suite are enabled, but if you do not need them, you can disable them during configuration like so:
```bash
cmake -S . -B build -DENABLE_TESTS=OFF -DENABLE_BENCHMARKS=OFF
```

Here is a brief summary of the CMake configuration options which can be enabled or disabled depending on the use case:
| Name              | Default Value | Description                                                                                                                                                        |
|-------------------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ENABLE_TESTS      | ON            | Enables the test suite.                                                                                                                                            |
| ENABLE_BENCHMARKS | ON            | Enables the benchmarks.                                                                                                                                            |
| USE_MPI           | OFF           | Enables the algorithms implemented with MPI (currently only Genetic Algorithm). Assumes that the variable CMAKE_CXX_COMPILER points to a MPI-capable C++ compiler. |

Summary of command line arguments and their purpose:
```bash
Command-line arguments:
 -h,  --help          Prints this message and exits.
 -a,  --algorithm     Sets the minimization algorithm to be used.
                      Must be one of: swarm_search, genetic_omp, genetic_mpi.
                      Default: genetic_omp.
 -d,  --dimensions    Sets the number of dimensions.
                      Must be >0.
                      Default: 2.
 -n,  --num-points    Sets the number of particles or creatures, depending on the algorithm.
                      Must be >0.
                      Default: 100.
 -i,  --iterations    Sets the maximum number of iterations.
                      Must be >0.
                      Default: 100.
 -lb, --lower-bound   Sets the lower boundary of the simulation space.
                      Must be finite.
                      Default: -100.
 -ub, --upper-bound   Sets the upper boundary of the simulation space.
                      Must be finite.
                      Default: 100.
 -sr, --survival-rate Sets the survival rate for the genetic algorithm.
                      Must be between 0.0 and 1.0.
                      Default: 0.5.
 -mr, --mutation-rate Sets the mutation rate for the genetic algorithm.
                      Must be between 0.0 and 1.0.
                      Default: 0.2.
 -f,  --function      Sets the function to be minimized.
                      Must be one of: sphere, euclideandistance, rosenbrock, rastrigin.
                      Default: sphere.
 -s,  --seed          Sets the seed for the random number generator.
                      Default: 4036653654.
 -j,  --jobs          Sets the number of threads to be used.
                      Must be >0 and <= 8.
                      Default: 1.
 -v, --verbose        Enables output on the console.
                      Does not need any arguments.
 -q, --quiet          Disables output on the console.
                      Does not need any arguments.

You can use it like so:
  ./build/main -d 2 -n 100 -i 100 -f sphere
```

### Note
Currently, there is no benchmark setup for the Genetic Algorithm implemented with MPI due to initialization issues originating from Google Benchmark. The recommended way to run benchmarks with the MPI algorithm, right now, is to compile with the `Release` configuration like so:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON
```
Then compiling the `main` executable and executing it specifying the algorithm on the command line like so:
```bash
cmake --build build -- main
mpirun -n 4 ./build/main -a genetic_mpi
```

## How to use this library
Inside the `example` folder you can find a simple example with CMake which declares a custom objective function and executes the OpenMP Genetic Algorithm.
This example automatically downloads the library from GitHub.
