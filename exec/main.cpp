#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <memory>
#include <random>
#include <omp.h>

#if defined(USE_MPI) && USE_MPI == 1
#include <mpi.h>

#include "DistributedGeneticAlgorithm.hpp"
#endif	// USE_MPI

#include "Swarm.hpp"
#include "GeneticAlgorithm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

enum class minimization_algorithm {
	SWARM_SEARCH,
	GENETIC_OPENMP
#if defined(USE_MPI) && USE_MPI == 1
	,
	GENETIC_MPI
#endif	// USE_MPI
};

std::ostream& operator<<(std::ostream& os, const minimization_algorithm algo) {
	switch (algo) {
		case minimization_algorithm::SWARM_SEARCH:
			os << "swarm_search";
			break;
		case minimization_algorithm::GENETIC_OPENMP:
			os << "genetic_omp";
			break;

#if defined(USE_MPI) && USE_MPI == 1
		case minimization_algorithm::GENETIC_MPI:
			os << "genetic_mpi";
			break;
#endif	// USE_MPI

		default:
			os << "Unknown";
			break;
	}
	return os;
}

void print_point(const size_t dimensions, const std::vector<double>& x, const double f) {
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
	std::cout << "  f(";
	for (size_t i = 0; i < dimensions; i++) {
		std::cout << std::scientific << x.at(i);
		if (i < dimensions - 1) {
			std::cout << ", ";
		}
	}
	std::cout << ") = " << std::scientific << f << std::endl;
	std::cout.copyfmt(oldState);
}

void run_swarm(const size_t dimensions, const size_t num_particles, const size_t max_iterations, const size_t seed,
			   const double lower_bound, const double upper_bound, const std::unique_ptr<ObjectiveFunction>& func,
			   const size_t n_threads) {
	std::vector<Particle> swarmParticles;

	for (size_t i = 0; i < num_particles; i++) {
		swarmParticles.push_back(Particle(dimensions, lower_bound, upper_bound, seed + i));
	}

	const double c1 = 2.0;
	const double c2 = 2.0;
	const double w_max = 0.9;
	const double w_min = 0.4;
	const double w = w_max;

	Swarm swarm = Swarm(swarmParticles, lower_bound, upper_bound, c1, c2, w, seed, *func, n_threads);

	const double beginning = omp_get_wtime();

	for (size_t i = 0; i < max_iterations; i++) {
		swarm.updateInertia(max_iterations, w_min, w_max);
		swarm.updateParticles();
		swarm.findBestFitness();
		std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
		std::cout << "  Current minimum: " << std::endl;
		print_point(dimensions, swarm.bestGlobalPosition, swarm.minimum);
		std::cout << std::endl;
	}

	const double end = omp_get_wtime();

	std::cout << std::endl;
	std::cout << "Minimum found:" << std::endl;
	print_point(dimensions, swarm.bestGlobalPosition, swarm.minimum);
	std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
			  << std::endl;
	std::cout << std::endl;
}

void run_genetic_openmp(const size_t dimensions, const size_t num_creatures, const size_t max_iterations,
						const size_t seed, const double lower_bound, const double upper_bound,
						const double mutation_rate, const double survival_rate,
						const std::unique_ptr<ObjectiveFunction>& func, const size_t n_threads) {
	std::vector<Creature> creatures;

	{
		std::mt19937 rnd{seed};
		std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
		for (size_t i{0}; i < num_creatures; i++) {
			std::vector<double> tmp(dimensions);
			std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
			creatures.push_back(Creature(tmp));
		}
	}

	GeneticAlgorithm ga(creatures, lower_bound, upper_bound, mutation_rate, survival_rate, *func, n_threads);

	const double beginning = omp_get_wtime();

	// Initial evaluation
	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();

		std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
		std::cout << "  Current minimum: " << std::endl;
		print_point(dimensions, ga.bestCreature.position, ga.bestCreature.fitness);
		std::cout << std::endl;
	}

	const double end = omp_get_wtime();

	std::cout << std::endl;
	std::cout << "Minimum found:" << std::endl;
	print_point(dimensions, ga.bestCreature.position, ga.bestCreature.fitness);
	std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
			  << std::endl;
	std::cout << std::endl;
}

#if defined(USE_MPI) && USE_MPI == 1
void run_genetic_mpi(const size_t dimensions, const size_t num_creatures, const size_t max_iterations,
					 const size_t seed, const double lower_bound, const double upper_bound, const double mutation_rate,
					 const double survival_rate, const std::unique_ptr<ObjectiveFunction>& func) {
	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	const size_t local_size = num_creatures / static_cast<size_t>(world_size);

	std::vector<std::vector<double>> creature_positions;

	// Only the "root process" initializes the creatures
	if (world_rank == 0) {
		std::mt19937 rnd{seed};
		std::uniform_real_distribution<double> dist{lower_bound, upper_bound};
		for (size_t i{0}; i < num_creatures; i++) {
			std::vector<double> tmp(dimensions);
			std::generate(tmp.begin(), tmp.end(), [&dist, &rnd]() { return dist(rnd); });
			creature_positions.push_back(tmp);
		}
	}

	// Split the creatures data among processes
	if (world_rank == 0) {
		// The root process sends to everyone
		for (int i = 1; i < world_size; i++) {
			for (size_t j{i * local_size}; j < (i + 1) * local_size; j++) {
				MPI_Send(creature_positions[j].data(), dimensions, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
	} else {
		// The other processes only receive
		creature_positions.resize(local_size);
		for (size_t i{0}; i < local_size; i++) {
			creature_positions[i].resize(dimensions);
			MPI_Recv(creature_positions[i].data(), dimensions, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	// The root process resizes its local vector
	if (world_rank == 0) {
		creature_positions.resize(local_size);
	}

	DistributedGeneticAlgorithm ga(world_rank, world_size, num_creatures, creature_positions, lower_bound, upper_bound,
								   mutation_rate, survival_rate, *func);

	MPI_Barrier(MPI_COMM_WORLD);
	const double beginning = MPI_Wtime();

	// Initial evaluation
	ga.evaluateCreatures();
	ga.sortCreatures();

	for (size_t i{0}; i < max_iterations; i++) {
		ga.applyCrossover(seed + i);
		ga.applyMutation(seed + i + 1);
		ga.evaluateCreatures();
		ga.sortCreatures();

		if (world_rank == 0) {
			std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
			std::cout << "  Current minimum: " << std::endl;
			print_point(dimensions, ga.creature_positions.at(ga.best_creature_index),
						ga.creature_fitnesses.at(ga.best_creature_index));
			std::cout << std::endl;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	const double end = MPI_Wtime();

	MPI_Finalize();

	if (world_rank == 0) {
		std::cout << std::endl;
		std::cout << "Minimum found:" << std::endl;
		print_point(dimensions, ga.creature_positions.at(ga.best_creature_index),
					ga.creature_fitnesses.at(ga.best_creature_index));
		std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
				  << std::endl;
		std::cout << std::endl;
	}
}
#endif	// USE_MPI

void die(const std::string& msg) {
	std::cerr << msg << std::endl;
	std::exit(-1);
}

size_t parse_integer(const std::string& arg, const std::string& short_arg_name, const std::string& long_arg_name) {
	size_t n;
	try {
		n = std::stoi(arg);
	} catch (const std::exception&) {
		die("Error: " + short_arg_name + ", " + long_arg_name + " requires a number.");
	}
	if (n <= 0) {
		die("Error: " + short_arg_name + ", " + long_arg_name + " must be >0.");
	}
	return n;
}

double parse_double(const std::string& arg, const std::string& short_arg_name, const std::string& long_arg_name) {
	double x;
	try {
		x = std::stod(arg);
	} catch (const std::exception&) {
		die("Error: " + short_arg_name + ", " + long_arg_name + " requires a number.");
	}
	if (!std::isfinite(x)) {
		die("Error: " + short_arg_name + ", " + long_arg_name + " must be finite.");
	}
	return x;
}

int main(const int argc, const char** argv) {
	std::random_device dev;

	minimization_algorithm algo = minimization_algorithm::GENETIC_OPENMP;
	size_t dimensions = 2;
	size_t num_points = 100;
	size_t max_iterations = 100;
	double lower_bound = -100.0;
	double upper_bound = 100.0;
	double mutation_rate = 0.2;
	double survival_rate = 0.5;
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
			std::cout << " -h,  --help          Prints this message and exits." << std::endl;
			std::cout << " -a,  --algorithm     Sets the minimization algorithm to be used." << std::endl;
			std::cout << "                      Must be one of: " << minimization_algorithm::SWARM_SEARCH << ", "
					  << minimization_algorithm::GENETIC_OPENMP
#if defined(USE_MPI) && USE_MPI == 1
					  << ", " << minimization_algorithm::GENETIC_MPI
#endif	// USE_MPI
					  << "." << std::endl;
			std::cout << "                      Default: " << algo << "." << std::endl;
			std::cout << " -d,  --dimensions    Sets the number of dimensions." << std::endl;
			std::cout << "                      Must be >0." << std::endl;
			std::cout << "                      Default: " << dimensions << "." << std::endl;
			std::cout << " -n,  --num-points    Sets the number of particles or creatures, depending on the algorithm."
					  << std::endl;
			std::cout << "                      Must be >0." << std::endl;
			std::cout << "                      Default: " << num_points << "." << std::endl;
			std::cout << " -i,  --iterations    Sets the maximum number of iterations." << std::endl;
			std::cout << "                      Must be >0." << std::endl;
			std::cout << "                      Default: " << max_iterations << "." << std::endl;
			std::cout << " -lb, --lower-bound   Sets the lower boundary of the simulation space." << std::endl;
			std::cout << "                      Must be finite." << std::endl;
			std::cout << "                      Default: " << lower_bound << "." << std::endl;
			std::cout << " -ub, --upper-bound   Sets the upper boundary of the simulation space." << std::endl;
			std::cout << "                      Must be finite." << std::endl;
			std::cout << "                      Default: " << upper_bound << "." << std::endl;
			std::cout << " -sr, --survival-rate Sets the survival rate for the genetic algorithm." << std::endl;
			std::cout << "                      Must be between 0.0 and 1.0." << std::endl;
			std::cout << "                      Default: " << survival_rate << "." << std::endl;
			std::cout << " -mr, --mutation-rate Sets the mutation rate for the genetic algorithm." << std::endl;
			std::cout << "                      Must be between 0.0 and 1.0." << std::endl;
			std::cout << "                      Default: " << mutation_rate << "." << std::endl;
			std::cout << " -f,  --function      Sets the function to be minimized." << std::endl;
			std::cout << "                      Must be one of: sphere, euclideandistance, rosenbrock, rastrigin."
					  << std::endl;
			std::cout << "                      Default: sphere." << std::endl;
			std::cout << " -s,  --seed          Sets the seed for the random number generator." << std::endl;
			std::cout << "                      Default: " << seed << "." << std::endl;
			std::cout << " -j,  --jobs          Sets the number of threads to be used." << std::endl;
			std::cout << "                      Must be >0 and <= " << omp_get_max_threads() << "." << std::endl;
			std::cout << "                      Default: " << n_threads << "." << std::endl;
			std::cout << std::endl;
			std::cout << "You can use it like so:" << std::endl;
			std::cout << "  " << argv[0] << " -d " << dimensions << " -n " << num_points << " -i " << max_iterations
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
			} else if (algo_name == "genetic_omp") {
				algo = minimization_algorithm::GENETIC_OPENMP;
			}
#if defined(USE_MPI) && USE_MPI == 1
			else if (algo_name == "genetic_mpi") {
				algo = minimization_algorithm::GENETIC_MPI;
			}
#endif	// USE_MPI
			else {
				die("Error: '" + algo_name + "' is not a known algorithm.");
			}
		} else if (arg == "-d" || arg == "--dimensions") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -d, --dimensions.");
			}
			dimensions = parse_integer(std::string(argv[i]), std::string("-d"), std::string("--dimensions"));
		} else if (arg == "-n" || arg == "--num-points") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -n, --num-points.");
			}
			num_points = parse_integer(std::string(argv[i]), std::string("-n"), std::string("--num-points"));
		} else if (arg == "-i" || arg == "--iterations") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -i, --iterations.");
			}
			max_iterations = parse_integer(std::string(argv[i]), std::string("-i"), std::string("--iterations"));
		} else if (arg == "-lb" || arg == "--lower-bound") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -lb, --lower-bound.");
			}
			lower_bound = parse_double(std::string(argv[i]), std::string("-lb"), std::string("--lower-bound"));
		} else if (arg == "-ub" || arg == "--upper-bound") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -ub, --upper-bound.");
			}
			upper_bound = parse_double(std::string(argv[i]), std::string("-ub"), std::string("--upper-bound"));
		} else if (arg == "-sr" || arg == "--survival-rate") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -sr, --survival-rate.");
			}
			survival_rate = parse_double(std::string(argv[i]), std::string("-sr"), std::string("--survival-rate"));
			if (survival_rate <= 0.0 || survival_rate >= 1.0) {
				die("Error: survival rate must be between 0.0 and 1.0.");
			}
		} else if (arg == "-mr" || arg == "--mutation-rate") {
			i++;
			if (i >= argc) {
				die("Error: missing argument for -mr, --mutation-rate.");
			}
			mutation_rate = parse_double(std::string(argv[i]), std::string("-mr"), std::string("--mutation-rate"));
			if (mutation_rate <= 0.0 || mutation_rate >= 1.0) {
				die("Error: mutation rate must be between 0.0 and 1.0.");
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
			n_threads = parse_integer(std::string(argv[i]), std::string("-j"), std::string("--jobs"));
			if (n_threads > static_cast<size_t>(omp_get_max_threads())) {
				die("Error: -j, --jobs must be <=" + std::to_string(omp_get_max_threads()) + ".");
			}
		} else {
			std::cerr << "Unknown argument '" << arg << "'" << std::endl;
			return -1;
		}
	}

	// Ensure that the lower bound is lower than the upper bound
	if (lower_bound > upper_bound) {
		std::swap(lower_bound, upper_bound);
	}

	if (algo == minimization_algorithm::SWARM_SEARCH) {
		run_swarm(dimensions, num_points, max_iterations, seed, lower_bound, upper_bound, func, n_threads);
	} else if (algo == minimization_algorithm::GENETIC_OPENMP) {
		run_genetic_openmp(dimensions, num_points, max_iterations, seed, lower_bound, upper_bound, mutation_rate,
						   survival_rate, func, n_threads);
	}
#if defined(USE_MPI) && USE_MPI == 1
	else if (algo == minimization_algorithm::GENETIC_MPI) {
		run_genetic_mpi(dimensions, num_points, max_iterations, seed, lower_bound, upper_bound, mutation_rate,
						survival_rate, func);
	}
#endif	// USE_MPI
	else {
		die("Error: unknown algorithm.");
	}

	return 0;
}
