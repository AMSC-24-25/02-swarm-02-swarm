#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <omp.h>

#if defined(USE_MPI) && USE_MPI == 1
#include <mpi.h>

#include "DistributedGeneticAlgorithm.hpp"
#include "DistributedDifferentialEvolution.hpp"
#include "DistributedSimulatedAnnealing.hpp"
#include "DistributedMultiStochasticTunnelling.hpp"
#endif	// USE_MPI

#include "Algorithm.hpp"

#include <Candidate.hpp>

#include "DifferentialEvolution.hpp"
#include <cassert>
#include "Particle.hpp"
#include "Swarm.hpp"
#include "GeneticAlgorithm.hpp"
#include "StochasticTunnelling.hpp"
#include "MultiStochasticTunnelling.hpp"
#include "Creature.hpp"
#include "Position.hpp"
#include "State.hpp"
#include "SimulatedAnnealing.hpp"
#include "FireflyAlgorithm.h"
#include "FireflyAlgorithm_Cuda.h"
#include "BFGSOptimizer.h"

namespace algorithm {

    namespace utils {

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

    }  // namespace utils

    std::pair<std::vector<double>, double> run_swarm(const size_t dimensions, const size_t num_particles,
                                                     const size_t max_iterations, const size_t seed,
                                                     const double lower_bound, const double upper_bound,
                                                     const std::unique_ptr<ObjectiveFunction>& func, const size_t n_threads,
                                                     const bool verbose) {
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

            if (verbose) {
                std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
                std::cout << "  Current minimum: " << std::endl;
                utils::print_point(dimensions, swarm.bestGlobalPosition, swarm.minimum);
                std::cout << std::endl;
            }
        }

        const double end = omp_get_wtime();

        if (verbose) {
            std::cout << std::endl;
            std::cout << "Minimum found:" << std::endl;
            utils::print_point(dimensions, swarm.bestGlobalPosition, swarm.minimum);
            std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
                      << std::endl;
            std::cout << std::endl;
        }

        return {swarm.bestGlobalPosition, swarm.minimum};
    };





   std::pair<std::vector<double>, double> run_multi_stochastic_tunnelling(const size_t dimensions,
                                                                           const size_t max_iterations, const size_t seed,
                                                                           const double lower_bound, const double upper_bound, const double sigma_max, const double sigma_min,
                                                                           const std::unique_ptr<ObjectiveFunction>& func, const double gamma,
                                                                           const double beta_adjust_factor, const bool verbose,double beta, const size_t tunnelling, const double beta_thresholding,
                                                                           const size_t num_positions, const size_t time_step_updating, const size_t n_threads) {

        std::vector<Position> pos;

        for(size_t i = 0; i < num_positions; i++){
            pos.push_back(Position(dimensions,lower_bound, upper_bound, seed, beta, *func, tunnelling, i));
        }


        MultiStochasticTunnelling stun = MultiStochasticTunnelling(pos, lower_bound, upper_bound, sigma_max, sigma_min, gamma, beta_adjust_factor, max_iterations, *func, beta_thresholding, num_positions, time_step_updating, dimensions, n_threads);


        const double beginning = omp_get_wtime();


        for(size_t i = 0; i < max_iterations; i++){
            stun.iteration(seed, i);
        }

        const double end = omp_get_wtime();

        if (verbose) {
            std::cout << std::endl;
            std::cout << "Minimum found:" << std::endl;
            utils::print_point(dimensions, stun.pos[0].best_position, stun.pos[0].f0);
            std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
                      << std::endl;
            std::cout << std::endl;
        }

        return {stun.pos[0].best_position, stun.pos[0].f0};


    }

    std::pair<std::vector<double>, double> run_differential_evolution(const size_t dimensions, const size_t num_candidates,
                                                                      const double lower_bound, const double upper_bound,
                                                                      const size_t seed, const size_t max_gen,
                                                                      const double F, const double CR,
                                                                      const std::unique_ptr<ObjectiveFunction>& func,
                                                                      const size_t n_threads, const bool verbose) {
        std::vector<Candidate> deCandidates;

        for (size_t i = 0; i < num_candidates; i++) {
            deCandidates.emplace_back(dimensions, lower_bound, upper_bound, seed + i,*func);
        }
        DifferentialEvolution de = DifferentialEvolution(deCandidates, dimensions,lower_bound,upper_bound,seed,max_gen,F,CR,*func,n_threads);
        de.findSol();
        const double beginning = omp_get_wtime();

        size_t iteration=0;
        while (iteration<max_gen) {
            iteration++;
            de.updateCandidate();
            de.findSol();
            const auto& best = de.getBestCandidate();
            if (verbose) {
                std::cout << "Iteration n. " << iteration << " / " << max_gen << std::endl;
                std::cout << "  Current minimum: " << std::endl;
                utils::print_point(dimensions, best.candidate, best.f0);
                std::cout << std::endl;
            }
        }

        const double end = omp_get_wtime();
        const auto& best = de.getBestCandidate();
        if (verbose) {
            std::cout << std::endl;
            std::cout << "Minimum found:" << std::endl;
            utils::print_point(dimensions, best.candidate, best.f0);
            std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
                      << std::endl;
            std::cout << std::endl;
        }

        return { best.candidate, best.f0 };
    }

std::pair<std::vector<double>, double> run_simulated_annealing(
                                                                   const size_t dimensions,
                                                                    const size_t max_iterations,
                                                                    const size_t dwell_iterations,
                                                                    const double initial_temperature,
                                                                    const double temperature_scale,
                                                                    const double initial_step_size,
                                                                    const double step_size_scale,
                                                                    const double boltzmann_constant,
                                                                    const std::vector<double>& initial_guess,
                                                                    const double lower_bound,
                                                                    const double upper_bound,
                                                                    const std::unique_ptr<ObjectiveFunction>& func,
                                                                    const size_t seed,
                                                                    const size_t n_threads,
                                                                    const bool verbose)
{
const size_t n_particles= n_threads; // Number of threads is equal to the number of particles

std::vector<std::unique_ptr<State>> bestStates(n_particles);
std::vector<double> bestCosts(n_threads, std::numeric_limits<double>::max());

const double start_time = omp_get_wtime();

#pragma omp parallel num_threads(n_threads)
{
    int tid = omp_get_thread_num();
    size_t local_seed = seed + tid;

    SimulatedAnnealing sa(
        *func,
        dimensions,
        max_iterations/n_particles,
        dwell_iterations/n_particles,
        initial_temperature,
        temperature_scale,
        initial_step_size,
        step_size_scale,
        boltzmann_constant,
        initial_guess,
        lower_bound,
        upper_bound,
        local_seed
    );

    sa.melt();
    sa.anneal();

    bestStates[tid] = std::make_unique<State>(sa.getBestState());
    bestCosts[tid] = bestStates[tid]->cost;
}
    const double end_time = omp_get_wtime();
    // Find minimum

    auto min_it = std::min_element(bestCosts.begin(), bestCosts.end());
    size_t best_idx = std::distance(bestCosts.begin(), min_it);

    

    if (verbose) {
        std::cout << "\nMinimum found:\n";
        algorithm::utils::print_point(dimensions, bestStates[best_idx]->values, bestStates[best_idx]->cost);
        std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end_time - start_time) << " seconds\n\n";
    }

    return {bestStates[best_idx]->values, bestCosts[best_idx]};

}

    std::pair<std::vector<double>, double> run_genetic_openmp(const size_t dimensions, const size_t num_creatures,
                                                              const size_t max_iterations, const size_t seed,
                                                              const double lower_bound, const double upper_bound,
                                                              const double mutation_rate, const double survival_rate,
                                                              const std::unique_ptr<ObjectiveFunction>& func,
                                                              const size_t n_threads, const bool verbose) {
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

            if (verbose) {
                std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
                std::cout << "  Current minimum: " << std::endl;
                utils::print_point(dimensions, ga.bestCreature.position, ga.bestCreature.fitness);
                std::cout << std::endl;
            }
        }

        const double end = omp_get_wtime();

        if (verbose) {
            std::cout << std::endl;
            std::cout << "Minimum found:" << std::endl;
            utils::print_point(dimensions, ga.bestCreature.position, ga.bestCreature.fitness);
            std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
                      << std::endl;
            std::cout << std::endl;
        }

        return {ga.bestCreature.position, ga.bestCreature.fitness};
    }






#ifdef USE_EIGEN
std::pair<std::vector<double>, double> run_firefly_bfgs(
	size_t dimensions,
	size_t num_fireflies,
	size_t max_iterations,
	size_t seed,
	double lower_bound,
	double upper_bound,
	const std::unique_ptr<ObjectiveFunction>& func,
	size_t n_threads,
	bool verbose,
	bool use_cuda,
	double alpha,
	double beta,
	double gamma
) {


	omp_set_num_threads(n_threads);

	auto objective = [&](const std::vector<double>& x) { return (*func)(x); };

	std::vector<double> best;
	double best_value = std::numeric_limits<double>::max();

	// 1 - FIRELY ALGORITHM
	if (use_cuda) {
#ifdef ENABLE_CUDA
		FireflyAlgorithm_Cuda algo(num_fireflies, dimensions, alpha, beta, gamma);
		algo.setObjectiveFunction(objective);
		best = algo.optimize(max_iterations);
		best_value = (*func)(best);
#else
		throw std::runtime_error("CUDA support not enabled in this build!");
#endif
	} else {
		FireflyAlgorithm algo(num_fireflies, dimensions, alpha, beta, gamma,lower_bound, upper_bound, seed);
		algo.setObjectiveFunction(objective);
		best = algo.optimize(max_iterations);
		best_value = (*func)(best);
	}

	// 2 - BFGS OPTIMIZER
	BFGSOptimizer bfgs(dimensions);
	bfgs.setObjective(objective);
	std::vector<double> refined = bfgs.optimize(best, max_iterations);
	double refined_value = (*func)(refined);

	if (verbose) {
		std::cout << "Firefly best: ";
		for (auto v : best) std::cout << v << " ";
		std::cout << "  value = " << best_value << std::endl;
		std::cout << "After BFGS: ";
		for (auto v : refined) std::cout << v << " ";
		std::cout << "  value = " << refined_value << std::endl;
	}

	// Restituisci la soluzione raffinata e il valore
	return {refined, refined_value};
}
#endif //USE_EIGEN


#if defined(USE_MPI) && USE_MPI == 1
    std::pair<std::vector<double>, double> run_genetic_mpi(const size_t dimensions, const size_t num_creatures,
                                                           const size_t max_iterations, const size_t seed,
                                                           const double lower_bound, const double upper_bound,
                                                           const double mutation_rate, const double survival_rate,
                                                           const std::unique_ptr<ObjectiveFunction>& func,
                                                           const bool verbose) {
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
                const size_t start = i * local_size;
                const size_t end = (i + 1) * local_size;
                for (size_t j{start}; j < end; j++) {
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

            if (verbose && world_rank == 0) {
                std::cout << "Iteration n. " << (i + 1) << " / " << max_iterations << std::endl;
                std::cout << "  Current minimum: " << std::endl;
                utils::print_point(dimensions, ga.creature_positions.at(ga.best_creature_index),
                                   ga.creature_fitnesses.at(ga.best_creature_index));
                std::cout << std::endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        const double end = MPI_Wtime();

        if (verbose && world_rank == 0) {
            std::cout << std::endl;
            std::cout << "Minimum found:" << std::endl;
            utils::print_point(dimensions, ga.creature_positions.at(ga.best_creature_index),
                               ga.creature_fitnesses.at(ga.best_creature_index));
            std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
                      << std::endl;
            std::cout << std::endl;
        }

        return {ga.creature_positions.at(ga.best_creature_index), ga.creature_fitnesses.at(ga.best_creature_index)};
    }


    std::pair<std::vector<double>, double> run_de_mpi(const size_t dimensions, const size_t num_candidates,
                                                      const double lower_bound, const double upper_bound,
                                                      const size_t seed, const size_t max_gen,
                                                      const double F, const double CR,
                                                      const std::unique_ptr<ObjectiveFunction>& func,
                                                      const bool verbose) {

        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        assert((int)num_candidates>=world_size);

        size_t base = num_candidates / world_size;
        size_t extras = num_candidates % world_size;
        // assegno ai primi processi gli extra, 1 extra ciascuno
        size_t my_size = base + (world_rank < (int)extras ? 1 : 0);

        std::vector<Candidate> deCandidates;
        deCandidates.reserve(my_size);

        size_t my_seed = seed + world_rank * 100000;
        // nb sto passando il seed non ho ancora creato la generazione!!!!!!!!
        // la generazione si crea dopo, se non faccio +i faccio generazioni con lo stesso motore e non cambia stato e genera sequenze uguali
        for (size_t i = 0; i < my_size; ++i) {
            deCandidates.emplace_back(dimensions, lower_bound, upper_bound, my_seed+i, *func);
        }

        DistributedDifferentialEvolution de(
                deCandidates, dimensions, lower_bound, upper_bound,
                my_seed, max_gen, F, CR, *func
        );


        std::vector<double> global_best_pos(dimensions);
        size_t iteration = 0;
        const double beginning = omp_get_wtime();

        std::pair<double,int> local_best;
        std::pair<double,int> global_best;

        while (iteration<max_gen) {
            iteration++;
            de.updateCandidate();
            de.findSol();
            local_best.first = de.bestCandidate->f0;
            local_best.second = world_rank;

            // otteniamo best fitness e rango del processo con best fitness
            // con queste due info e con la local_best_position di ogni processo
            // poi possiamo ricavare anche la global_best_position
            MPI_Allreduce(
                    &local_best,
                    &global_best,
                    1,
                    MPI_DOUBLE_INT,
                    MPI_MINLOC,
                    MPI_COMM_WORLD
            );

            if(world_rank == global_best.second) {
                std::copy(de.bestCandidate->candidate.begin(), de.bestCandidate->candidate.end(),
                          global_best_pos.begin());
            }
            // broadcastiamo a tutti i processi la best position
            MPI_Bcast(
                    global_best_pos.data(),
                    (int)dimensions,
                    MPI_DOUBLE,
                    global_best.second,
                    MPI_COMM_WORLD
            );

            if (verbose && world_rank == global_best.second) {
                std::cout << "Iteration n. " << iteration << " / " << max_gen << std::endl;
                std::cout << "  Current minimum: " << std::endl;
                utils::print_point(dimensions, de.bestCandidate->candidate, de.bestCandidate->f0);
                std::cout << std::endl;
            }


        }

        const double end = omp_get_wtime();

        if (verbose && world_rank==0) {
            std::cout << std::endl;
            std::cout << "Minimum found:" << std::endl;
            utils::print_point(dimensions, global_best_pos , global_best.first);
            std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
                      << std::endl;
            std::cout << std::endl;
        }

        return {global_best_pos, global_best.first};
    }

std::pair<std::vector<double>, double> run_sa_mpi(const size_t dimensions,
												const size_t max_iterations, const size_t dwell_iterations,
												const double initial_temperature, const double temperature_scale,
												const double initial_step_size, const double step_size_scale,
												const double boltzmann_constant,
												const double lower_bound,const double upper_bound,
												const std::unique_ptr<ObjectiveFunction>& func, const size_t seed,
											    const bool verbose) {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   // Inizializzazione randomica diversa per ogni processo
    std::vector<double> my_guess(dimensions);
    std::mt19937 local_gen(seed + world_rank * 99991);
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    for (size_t i = 0; i < dimensions; ++i) {
        my_guess[i] = dist(local_gen);
    }

    // Costruzione oggetto SA locale
    DistributedSimulatedAnnealing sa(*func,
                                     dimensions,
                                     max_iterations,
                                     dwell_iterations,
                                     initial_temperature,
                                     temperature_scale,
                                     initial_step_size,
                                     step_size_scale,
                                     boltzmann_constant,
                                     my_guess,
                                     lower_bound,
                                     upper_bound,
                                     seed + world_rank * 1000);

    const double start_time = omp_get_wtime();

    sa.melt();    // fase di melt
    sa.anneal();  // fase di raffreddamento

    const State& local_best = sa.getBestState();
    std::pair<double, int> local_result = {local_best.cost, world_rank};
    std::pair<double, int> global_result;

    // Trova la migliore soluzione globale
    MPI_Allreduce(
        &local_result,
        &global_result,
        1,
        MPI_DOUBLE_INT,
        MPI_MINLOC,
        MPI_COMM_WORLD
    );

    std::vector<double> global_best_position(dimensions);
    if (world_rank == global_result.second) {
        std::copy(local_best.values.begin(), local_best.values.end(), global_best_position.begin());
    }

    // Trasmette la miglior soluzione a tutti
    MPI_Bcast(global_best_position.data(), static_cast<int>(dimensions), MPI_DOUBLE, global_result.second, MPI_COMM_WORLD);

    const double end_time = omp_get_wtime();

    if (verbose && world_rank == 0) {
        std::cout << "Best solution found by process " << global_result.second << ":" << std::endl;
        utils::print_point(dimensions, global_best_position, global_result.first);
        std::cout << "Execution time: " << std::fixed << std::setprecision(6) << (end_time - start_time) << " seconds" << std::endl;
    }

    return {global_best_position, global_result.first};
}





    std::pair<std::vector<double>, double> run_tunnelling_mpi(const size_t dimensions,
                                                              const size_t max_iterations, const size_t seed,
                                                              const double lower_bound, const double upper_bound, const double sigma_max, const double sigma_min,
                                                              const  ObjectiveFunction& func, const double gamma,
                                                              const double beta_adjust_factor, const bool verbose,double beta, const size_t tunnelling, const double beta_thresholding,
                                                              const size_t num_positions, const size_t time_step_updating) {
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);




        std::vector<Position> pos;

        if (world_rank == 0) {
            for(size_t i = 0; i < num_positions; i++){
                pos.push_back(Position(dimensions,lower_bound, upper_bound, seed, beta, func, tunnelling, i));
            }
        }


        DistributedMultiStochasticTunnelling stun = DistributedMultiStochasticTunnelling(world_rank, world_size, pos, lower_bound, upper_bound, sigma_max, sigma_min, gamma, beta_adjust_factor, max_iterations, func, beta_thresholding, num_positions, time_step_updating, dimensions);




        MPI_Barrier(MPI_COMM_WORLD);
        const double beginning = MPI_Wtime();

        for(size_t i = 0; i < max_iterations; i++){
            if(i == floor(max_iterations*(1.0/3.0))){
                stun.update_beta_thresholding(1);
            }else if(i == floor(max_iterations*(2.0/3.0))){
                stun.update_beta_thresholding(2);
            }
            stun.iteration(seed, i);

        }

        MPI_Barrier(MPI_COMM_WORLD);
        const double end = MPI_Wtime();

        if (verbose && world_rank == 0) {
            std::cout << std::endl;
            std::cout << "Minimum found:" << std::endl;
            utils::print_point(dimensions, stun.pos[0].best_position, stun.pos[0].f0);
            std::cout << "  Total execution time: " << std::fixed << std::setprecision(6) << (end - beginning) << " seconds"
                      << std::endl;
            std::cout << std::endl;
        }

        return {stun.pos[0].best_position, stun.pos[0].f0};

    }

#endif	// USE_MPI

}  // namespace algorithm

