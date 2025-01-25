#ifndef DISTRIBUTED_GENETIC_ALGORITHM_HPP
#define DISTRIBUTED_GENETIC_ALGORITHM_HPP

#include <vector>

#include "Creature.hpp"
#include "ObjectiveFunction.hpp"

class DistributedGeneticAlgorithm {
   public:
	std::vector<Creature> creatures;
	const double lower_bound;
	const double upper_bound;
	const double survival_rate;
	const double mutation_rate;
	Creature bestCreature;
	const ObjectiveFunction& func;

	DistributedGeneticAlgorithm(const std::vector<Creature>& creatures, const double lower_bound,
								const double upper_bound, const double mutation_rate, const double survival_rate,
								ObjectiveFunction& func);

	void evaluateCreatures();

	void sortCreatures();

	void applyCrossover(const size_t seed);

	void applyMutation(const size_t seed);
};

#endif	// DISTRIBUTED_GENETIC_ALGORITHM_HPP
