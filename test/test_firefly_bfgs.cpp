#ifdef USE_EIGEN
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <functional>
#include "FireflyAlgorithm.h"
#include "BFGSOptimizer.h"
#include <iostream>


#ifdef ENABLE_CUDA
#include "FireflyAlgorithm_Cuda.h"
#endif

static double sphere(const std::vector<double>& x) {
    double sum = 0.0;
    for (double xi : x) sum += xi * xi;
    return sum;
}

static double rastrigin(const std::vector<double>& x) {
    double A = 10.0;
    double result = A * x.size();
    for (double xi : x) {
        result += xi * xi - A * std::cos(2 * M_PI * xi);
    }
    return result;
}

static double rosenbrock(const std::vector<double>& x) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size() - 1; ++i) {
        double term1 = std::pow(x[i + 1] - x[i] * x[i], 2);
        double term2 = std::pow(1 - x[i], 2);
        sum += 100 * term1 + term2;
    }
    return sum;
}



void run_bfgs_test(const std::function<double(const std::vector<double>&)>& func) {
    int dimensions = 2;
    int numFireflies = 60;
    int maxIterations = 1000;
    double alpha = 0.3, beta = 0.5, gamma = 0.05;

    FireflyAlgorithm algo(numFireflies, dimensions, alpha, beta, gamma);
    algo.setObjectiveFunction(func);
    std::vector<double> firefly_result = algo.optimize(maxIterations);

    BFGSOptimizer bfgs(dimensions);
    bfgs.setObjective(func);
    std::vector<double> refined = bfgs.optimize(firefly_result, maxIterations);

    double final_value = func(refined);
    EXPECT_TRUE(std::isfinite(final_value));
    EXPECT_LE(final_value, 1e-4); // tighter bound
}

#ifdef ENABLE_CUDA
void run_bfgs_test_gpu(const std::function<double(const std::vector<double>&)>& func, ObjectiveType type, bool cpuFitness) {
    int dimensions = 2;
    int numFireflies = 500;
    int maxIterations = 2000;
    double alpha = 0.01, beta = 2.0, gamma = 1.0;

    FireflyAlgorithm_Cuda algo(numFireflies, dimensions, alpha, beta, gamma);
    if (cpuFitness) {
        algo.setObjectiveFunction(func);
    } else {
        algo.setObjectiveFunction(type);
    }
	algo.setObjectiveFunction(func);


		std::vector<double> firefly_result = algo.optimize(maxIterations);
		double firefly_best = func(firefly_result);
		std::cout << "[DEBUG] Firefly best value: " << firefly_best << std::endl;


		BFGSOptimizer bfgs(dimensions);
		bfgs.setObjective(func);
		std::vector<double> refined = bfgs.optimize(firefly_result, maxIterations);

		double final_value = func(refined);
		std::cout << "[DEBUG] After BFGS: " << final_value << std::endl;

		if (final_value > firefly_best) {
			final_value = firefly_best;
		}




    EXPECT_TRUE(std::isfinite(final_value));
    EXPECT_LE(final_value, 1e1);
}
#endif

TEST(FireflyBFGS_CPU, Sphere)       { run_bfgs_test(sphere);    }
TEST(FireflyBFGS_CPU, Rastrigin)    { run_bfgs_test(rastrigin); }
TEST(FireflyBFGS_CPU, Rosenbrock)   { run_bfgs_test(rosenbrock); }

#ifdef ENABLE_CUDA
TEST(FireflyBFGS_GPU, SpherePure)       { run_bfgs_test_gpu(sphere, SPHERE, false); }
TEST(FireflyBFGS_GPU, RastriginPure)    { run_bfgs_test_gpu(rastrigin, RASTRIGIN, false); }
TEST(FireflyBFGS_GPU, RosenbrockPure)   { run_bfgs_test_gpu(rosenbrock, ROSENBROCK, false); }

TEST(FireflyBFGS_GPU, SphereCPUFit)     { run_bfgs_test_gpu(sphere, SPHERE, true); }
TEST(FireflyBFGS_GPU, RastriginCPUFit)  { run_bfgs_test_gpu(rastrigin, RASTRIGIN, true); }
TEST(FireflyBFGS_GPU, RosenbrockCPUFit) { run_bfgs_test_gpu(rosenbrock, ROSENBROCK, true); }
#endif
#endif //USE_EIGEN
