#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <vector>

#include "Creature.hpp"
#include "ObjectiveFunction.hpp"

class GeneticAlgorithm {
   private:
	const double lower_bound;
	const double upper_bound;
	const double survival_rate;
	const double mutation_rate;
	const size_t n_threads;
	const ObjectiveFunction& func;

   public:
	std::vector<Creature> creatures;
	Creature bestCreature;

	GeneticAlgorithm(const std::vector<Creature>& creatures, const double lower_bound, const double upper_bound,
					 const double mutation_rate, const double survival_rate, ObjectiveFunction& func,
					 const size_t n_threads);

	void evaluateCreatures();

	void sortCreatures();

	void applyCrossover(const size_t seed);

	void applyMutation(const size_t seed);
};

#endif	// GENETIC_ALGORITHM_HPP
