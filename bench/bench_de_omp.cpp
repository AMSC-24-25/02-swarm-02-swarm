#include <benchmark/benchmark.h>

#include <random>
#include <chrono>

#include "Algorithm.hpp"
#include "Rosenbrock.hpp"

static void DeOpenMP_Rosenbrock(benchmark::State& state) {
    const size_t dimensions = 2;
    const size_t num_creatures = state.range(0);
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
                                                      F, CR, r, 1, false);

        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();

        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.SetIterationTime(elapsed_seconds.count());
        state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
    }
}

BENCHMARK(DeOpenMP_Rosenbrock)->Range(4, 1'000)->RangeMultiplier(2)->UseManualTime();
