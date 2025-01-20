#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <memory>
#include <random>
#include <omp.h>

#include "Swarm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

enum class minimization_algorithm { SWARM_SEARCH, GENETIC };

std::ostream& operator<<(std::ostream& os, const minimization_algorithm algo) {
	switch (algo) {
		case minimization_algorithm::SWARM_SEARCH:
			os << "swarm_search";
			break;
		case minimization_algorithm::GENETIC:
			os << "genetic";
			break;
		default:
			os << "Unknown";
			break;
	}
	return os;
}

void print_minimum(const Swarm& swarm, const size_t dimensions) {
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
	std::cout << "  f(";
	for (size_t i = 0; i < dimensions; i++) {
		std::cout << std::scientific << swarm.bestGlobalPosition[i];
		if (i < dimensions - 1) {
			std::cout << ", ";
		}
	}
	std::cout << ") = " << std::scientific << swarm.minimum << std::endl;
	std::cout.copyfmt(oldState);
}

void run_swarm(const int dimensions, const int num_particles, const int max_iterations, const size_t seed,
			   const std::unique_ptr<ObjectiveFunction>& func, const size_t n_threads) {
	std::vector<Particle> swarmParticles;

	std::vector<double> lowerBound(dimensions, -100.0);
	std::vector<double> upperBound(dimensions, 100.0);

	for (int i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lowerBound, upperBound, seed));
	}

	const double c1 = 2.0;
	const double c2 = 2.0;
	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Swarm swarm = Swarm(swarmParticles, lowerBound, upperBound, c1, c2, w, seed, *func, n_threads);

	const auto beginning = omp_get_wtime();

	for (int i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();
		std::cout << "Iteration n. " << i << " / " << max_iterations << std::endl;
		std::cout << "  Current minimum: " << std::endl;
		print_minimum(swarm, dimensions);
		std::cout << std::endl;
	}

	const auto end = omp_get_wtime();

	std::cout << std::endl;
	std::cout << "Minimum found:" << std::endl;
	print_minimum(swarm, dimensions);
	std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
			  << std::endl;
	std::cout << std::endl;
}

void run_genetic() {
	std::cout << "TODO: implement genetic algorithm." << std::endl;
}

void die(const std::string& msg) {
	std::cerr << msg << std::endl;
	std::exit(-1);
}

int main(const int argc, const char** argv) {
	std::random_device dev;

	minimization_algorithm algo = minimization_algorithm::GENETIC;
	int dimensions = 2;
	int num_particles = 100;
	int max_iterations = 100;
	std::unique_ptr<ObjectiveFunction> func = std::make_unique<Sphere>();
	size_t seed = dev();
	size_t n_threads = 1;

	for (int i = 1; i < argc; i++) {
		std::string arg = std::string(argv[i]);
		if (arg == "-h" || arg == "--help") {
			std::cout << std::endl;
			std::cout
				<< "  Welcome to group 02 project for Advanced Methods for Scientific Computing course a.y. 2024/2025"
				<< std::endl;
			std::cout << "       SWARM SEARCH" << std::endl;
			std::cout << std::endl;
			std::cout << "Command-line arguments:" << std::endl;
			std::cout << " -h, --help        Prints this message and exits." << std::endl;
			std::cout << " -a, --algorithm   Sets the minimization algorithm to be used. Must be one of: "
					  << minimization_algorithm::SWARM_SEARCH << ", " << minimization_algorithm::GENETIC
					  << ". Default: " << algo << "." << std::endl;
			std::cout << " -d, --dimensions  Sets the number of dimensions. Must be >0. Default: " << dimensions << "."
					  << std::endl;
			std::cout << " -p, --particles   Sets the number of particles. Must be >0. Default: " << num_particles
					  << "." << std::endl;
			std::cout << " -i, --iterations  Sets the number of iterations. Must be >0. Default: " << max_iterations
					  << "." << std::endl;
			std::cout << " -f, --function    Sets the function to be minimized. Must be one of: sphere, "
						 "euclideandistance, rosenbrock, rastrigin. Default: sphere."
					  << std::endl;
			std::cout << " -s, --seed        Sets the seed for the random number generator. Default: " << seed << "."
					  << std::endl;
			std::cout << " -j, --jobs        Sets the number of threads to be used. Must be >0 and <= "
					  << omp_get_max_threads() << ". Default: " << n_threads << "." << std::endl;
			std::cout << std::endl;
			std::cout << "You can use it like so:" << std::endl;
			std::cout << "  " << argv[0] << " -d " << dimensions << " -p " << num_particles << " -i " << max_iterations
					  << " -f sphere" << std::endl;
			std::cout << std::endl;
			return 0;
		} else if (arg == "-a" || arg == "--algorithm") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -a, --algorithm.");
			}
			std::string algo_name = std::string(argv[i]);
			std::transform(algo_name.begin(), algo_name.end(), algo_name.begin(),
						   [](const char ch) { return std::tolower(ch); });
			if (algo_name == "swarm_search") {
				algo = minimization_algorithm::SWARM_SEARCH;
			} else if (algo_name == "genetic") {
				algo = minimization_algorithm::GENETIC;
			} else {
				die("Error: '" + algo_name + "' is not a known algorithm.");
			}
		} else if (arg == "-d" || arg == "--dimensions") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -d, --dimensions.");
			}
			try {
				dimensions = std::stoi(std::string(argv[i]));
			} catch (const std::exception&) {
				die("Error: -d, --dimensions requires a number.");
			}
			if (dimensions <= 0) {
				die("Error: -d, --dimensions must be >0.");
			}
		} else if (arg == "-p" || arg == "--particles") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -p, --particles.");
			}
			try {
				num_particles = std::stoi(std::string(argv[i]));
			} catch (const std::exception&) {
				die("Error: -p, --particles requires a number.");
			}
			if (num_particles <= 0) {
				die("Error: -p, --particles must be >0.");
			}
		} else if (arg == "-i" || arg == "--iterations") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -i, --iterations.");
			}
			try {
				max_iterations = std::stoi(std::string(argv[i]));
			} catch (const std::exception&) {
				die("Error: -i, --iterations requires a number.");
			}
			if (max_iterations <= 0) {
				die("Error: -i, --iterations must be >0.");
			}
		} else if (arg == "-f" || arg == "--function") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -f, --function.");
			}
			std::string function_name = std::string(argv[i]);
			std::transform(function_name.begin(), function_name.end(), function_name.begin(),
						   [](const char ch) { return std::tolower(ch); });
			if (function_name == "sphere") {
				func = std::make_unique<Sphere>();
			} else if (function_name == "euclideandistance") {
				func = std::make_unique<EuclideanDistance>();
			} else if (function_name == "rosenbrock") {
				func = std::make_unique<Rosenbrock>();
			} else if (function_name == "rastrigin") {
				func = std::make_unique<Rastrigin>();
			} else {
				die("Error: '" + function_name + "' is not a known function.");
			}
		} else if (arg == "-s" || arg == "--seed") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -s, --seed.");
			}
			try {
				seed = std::stoi(std::string(argv[i]));
			} catch (const std::exception&) {
				die("Error: -s, --seed requires a number.");
			}
		} else if (arg == "-j" || arg == "--jobs") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -j, --jobs.");
			}
			try {
				n_threads = std::stoi(std::string(argv[i]));
			} catch (const std::exception&) {
				die("Error: -j, --jobs requires a number.");
			}
			if (n_threads <= 0 || n_threads > static_cast<size_t>(omp_get_max_threads())) {
				die("Error: -j, --jobs must be >0 and <=" + std::to_string(omp_get_max_threads()) + ".");
			}
		} else {
			std::cerr << "Unknown argument '" << arg << "'" << std::endl;
			return -1;
		}
	}

	if (algo == minimization_algorithm::SWARM_SEARCH) {
		run_swarm(dimensions, num_particles, max_iterations, seed, func, n_threads);
	} else if (algo == minimization_algorithm::GENETIC) {
		run_genetic();
	} else {
		die("Error: unknown algorithm.");
	}

	return 0;
}
