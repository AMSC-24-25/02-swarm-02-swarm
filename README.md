# Project for AMSC 2024/2025 course
This is the project for the group 02. It implements some algorithms:
 - [Swarm Search](https://github.com/AMSC-24-25/02-swarm-02-swarm/blob/main/SwarmSearch.md)
 - [Genetic Algorithm](https://github.com/AMSC-24-25/02-swarm-02-swarm/blob/main/GeneticAlgorithm.md)

## How to compile and run
This project uses CMake as build system.
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

By default, the test suite and the benchmark suite are enabled, but if you do not need them, you can disable them during configuration like so:
```bash
cmake -S . -B build -DENABLE_TESTS=OFF -DENABLE_BENCHMARKS=OFF
```
