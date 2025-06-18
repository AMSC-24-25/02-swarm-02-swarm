#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <omp.h>

#if defined(USE_MPI) && USE_MPI == 1
#include <mpi.h>

#endif	// USE_MPI

#include "Algorithm.hpp"
#include "Sphere.hpp"
#include "EuclideanDistance.hpp"
#include "Rosenbrock.hpp"
#include "Rastrigin.hpp"

enum class minimization_algorithm {
	SWARM_SEARCH,
	GENETIC_OPENMP,
    DE_OPENMP,
	SA_OPENMP,

#if defined(USE_MPI) && USE_MPI == 1
	GENETIC_MPI,
    DE_MPI,
	SA_MPI
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

        case minimization_algorithm::DE_OPENMP:
            os << "differential_omp";
            break;
		
		case minimization_algorithm::SA_OPENMP:
            os << "simulated_omp";
            break;	

#if defined(USE_MPI) && USE_MPI == 1
		case minimization_algorithm::GENETIC_MPI:
			os << "genetic_mpi";
			break;
        case minimization_algorithm::DE_MPI:
            os << "differential_mpi";
            break;
		case minimization_algorithm::SA_MPI:
            os << "simulated_mpi";
            break;
#endif	// USE_MPI

		default:
			os << "Unknown";
			break;
	}
	return os;
}

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

    /* -------------------------------*/
    /* GENETIC PARAMETERS*/
	double mutation_rate = 0.2;
	double survival_rate = 0.5;
    /* DIFFERENTIAL EVOLUTION PARAMETERS*/
    double F = 0.5;
    double CR = 0.8;
	/* SIMULATED ANNEALING PARAMETERS*/
    int dwell=1000;
	double initial_temperature = 10.0;
	double temperature_scale = 0.98;
	double initial_step_size = 1.0;
	double step_size_scale = 0.98;
	double boltzmann_constant = 1.0;
	std::vector<double> initial_guess(dimensions, 5.0);
    /* -------------------------------*/

	std::unique_ptr<ObjectiveFunction> func = std::make_unique<Sphere>();
	size_t seed = dev();
	size_t n_threads = 1;
	bool verbose = true;

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
					  << minimization_algorithm::GENETIC_OPENMP << ", " << minimization_algorithm::DE_OPENMP << ", " 
					  << minimization_algorithm::SA_OPENMP
#if defined(USE_MPI) && USE_MPI == 1
					  << ", " << minimization_algorithm::GENETIC_MPI << ", " << minimization_algorithm::DE_MPI  << ", "
	   << minimization_algorithm::SA_MPI
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
			std::cout << " -v, --verbose        Enables output on the console." << std::endl;
			std::cout << "                      Does not need any arguments." << std::endl;
			std::cout << " -q, --quiet          Disables output on the console." << std::endl;
			std::cout << "                      Does not need any arguments." << std::endl;
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
            } else if (algo_name == "differential_omp"){
                algo = minimization_algorithm::DE_OPENMP;
            }else if (algo_name == "simulated_omp") {
				algo = minimization_algorithm::SA_OPENMP;
			}
#if defined(USE_MPI) && USE_MPI == 1
			else if (algo_name == "genetic_mpi") {
				algo = minimization_algorithm::GENETIC_MPI;
			}else if(algo_name == "differential_mpi"){
                algo = minimization_algorithm::DE_MPI;
            }else if (algo_name == "simulated_mpi") {
	            algo = minimization_algorithm::SA_MPI;
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
		} else if (arg == "-v" || arg == "--verbose") {
			verbose = true;
		} else if (arg == "-q" || arg == "--quiet") {
			verbose = false;
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
		algorithm::run_swarm(dimensions, num_points, max_iterations, seed, lower_bound, upper_bound, func, n_threads,
							 verbose);
	} else if (algo == minimization_algorithm::GENETIC_OPENMP) {
		algorithm::run_genetic_openmp(dimensions, num_points, max_iterations, seed, lower_bound, upper_bound,
									  mutation_rate, survival_rate, func, n_threads, verbose);
	} else if(algo == minimization_algorithm::DE_OPENMP){
        algorithm::run_differential_evolution(dimensions,num_points,lower_bound,upper_bound,seed,max_iterations,F,CR,func,n_threads,verbose);
    }else if (algo == minimization_algorithm::SA_OPENMP) {
  algorithm::run_simulated_annealing(dimensions, max_iterations, dwell, initial_temperature,
									 temperature_scale, initial_step_size, step_size_scale, boltzmann_constant,	initial_guess,
									lower_bound, upper_bound, func, seed, n_threads, verbose);}
#if defined(USE_MPI) && USE_MPI == 1
	else if (algo == minimization_algorithm::GENETIC_MPI) {
		MPI_Init(NULL, NULL);
		algorithm::run_genetic_mpi(dimensions, num_points, max_iterations, seed, lower_bound, upper_bound,
								   mutation_rate, survival_rate, func, verbose);
		MPI_Finalize();
	}else if (algo == minimization_algorithm::DE_MPI){
        MPI_Init(NULL,NULL);
        algorithm::run_de_mpi(dimensions,num_points,lower_bound,upper_bound,seed,max_iterations,F,CR,func,verbose);
        MPI_Finalize();
    }else if( algo == minimization_algorithm::SA_MPI) {
		MPI_Init(NULL, NULL);
  		algorithm::run_sa_mpi(dimensions, max_iterations, dwell, initial_temperature, temperature_scale,
		initial_step_size, step_size_scale, boltzmann_constant, lower_bound, upper_bound, func, seed, verbose);
  		MPI_Finalize();
   }
#endif	// USE_MPI
	else {
		die("Error: unknown algorithm.");
	}

	return 0;
}
