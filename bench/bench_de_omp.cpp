#include <benchmark/benchmark.h>

#include <random>
#include <chrono>
#include "Algorithm.hpp"
#include "Rosenbrock.hpp"
#include "omp.h"

static void DeOpenMP_Rosenbrock(benchmark::State& state) {
    const size_t dimensions = 5;
    const size_t num_creatures = state.range(0);
    const int    num_threads   = state.threads();   // prendo i thread da ThreadRange
    omp_set_num_threads(num_threads);

    const size_t max_iterations = 100;
    const size_t seed = 42;
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const double F = 0.5;
    const double CR = 0.5;
    const std::unique_ptr<ObjectiveFunction> r = std::make_unique<Rosenbrock>();

    for (auto _ : state) {
        const auto start = std::chrono::high_resolution_clock::now();

        std::pair<std::vector<double>, double> result =
                algorithm::run_differential_evolution(dimensions, num_creatures, lower_bound, upper_bound,seed,max_iterations,
                                                      F, CR, r, num_threads, false);

        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();

        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.SetIterationTime(elapsed_seconds.count());
        state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
    }
}

BENCHMARK(DeOpenMP_Rosenbrock)
    ->ArgName("creatures")         // usa un solo argomento di range
    ->RangeMultiplier(2)
    ->Range(4, 1000000)             // creatures = 4, 8, 16, …, 65536
    ->ThreadRange(1, 16)           // threads   = 1, 2, 4, 8, 16
    ->UseManualTime();

