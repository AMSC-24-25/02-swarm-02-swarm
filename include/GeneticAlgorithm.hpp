#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <vector>

#include "Creature.hpp"
#include "ObjectiveFunction.hpp"

class GeneticAlgorithm {
   public:
	std::vector<Creature> creatures;
	const double lower_bound;
	const double upper_bound;
	const double survival_rate;
	const double mutation_rate;
	Creature bestCreature;
	const ObjectiveFunction& func;

	GeneticAlgorithm(const std::vector<Creature>& creatures, const double lower_bound, const double upper_bound,
					 const double mutation_rate, const double survival_rate, ObjectiveFunction& func);

	void evaluateCreatures();

	void sortCreatures();

	void applyCrossover();

	void applyMutation();
};

#endif	// GENETIC_ALGORITHM_HPP
