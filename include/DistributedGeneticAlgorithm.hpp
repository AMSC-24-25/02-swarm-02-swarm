#ifndef DISTRIBUTED_GENETIC_ALGORITHM_HPP
#define DISTRIBUTED_GENETIC_ALGORITHM_HPP

#include <vector>
#include <cstddef>

#include "ObjectiveFunction.hpp"

class DistributedGeneticAlgorithm {
   private:
	const int world_rank;
	const int world_size;
	const size_t total_creatures;
	const double lower_bound;
	const double upper_bound;
	const double survival_rate;
	const double mutation_rate;
	const ObjectiveFunction& func;

   public:
	std::vector<std::vector<double>> creature_positions;
	std::vector<double> creature_fitnesses;
	size_t best_creature_index;

	DistributedGeneticAlgorithm(const int world_rank, const int world_size, const size_t total_creatures,
								const std::vector<std::vector<double>>& creature_positions, const double lower_bound,
								const double upper_bound, const double mutation_rate, const double survival_rate,
								ObjectiveFunction& func);

	void evaluateCreatures();

	void sortCreatures();

	void applyCrossover(const size_t seed);

	void applyMutation(const size_t seed);
};

#endif	// DISTRIBUTED_GENETIC_ALGORITHM_HPP
