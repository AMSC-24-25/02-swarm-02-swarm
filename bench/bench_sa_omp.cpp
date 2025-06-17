#include <benchmark/benchmark.h>
#include <chrono>
#include <vector>
#include <memory>

#include "Algorithm.hpp"
#include "Rosenbrock.hpp"

static void SAOpenMP_Rosenbrock(benchmark::State& state) {
    const size_t dimensions = 2;

    const size_t n_threads = state.range(0);
    const size_t max_iterations = state.range(1);

    const size_t dwell_iterations = 20;
    const double initial_temperature = 100.0;
    const double temperature_scale = 0.9;
    const double initial_step_size = 0.5;
    const double step_size_scale = 0.95;
    const double boltzmann_constant = 1.0;

    const std::vector<double> initial_guess(dimensions, 0.0);  // esempio: (0.0, 0.0)
    const double lower_bound = -10.0;
    const double upper_bound = 10.0;
    const size_t seed = 42;

    const std::unique_ptr<ObjectiveFunction> func = std::make_unique<Rosenbrock>();

    for (auto _ : state) {
        const auto start = std::chrono::high_resolution_clock::now();

        auto result = algorithm::run_simulated_annealing(
            dimensions,
            max_iterations,
            dwell_iterations,
            initial_temperature,
            temperature_scale,
            initial_step_size,
            step_size_scale,
            boltzmann_constant,
            initial_guess,
            lower_bound,
            upper_bound,
            func,
            seed,
            n_threads,
            false // verbose
        );

        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();

        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.SetIterationTime(elapsed_seconds.count());
        state.counters["Iteration/s"] = static_cast<double>(max_iterations) / elapsed_seconds.count();
        state.counters["BestCost"] = result.second;
    }
}

BENCHMARK(SAOpenMP_Rosenbrock)
    ->ArgNames({"threads", "max_iterations"})
    ->Ranges({{1, 8},         // threads: 1, 2, 4, 8
              {100, 800}})    // iterations: 100, 200, 400, 800
    ->RangeMultiplier(2)
    ->UseManualTime();
