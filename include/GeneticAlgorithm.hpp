#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <vector>

class GeneticAlgorithm {
   public:
	using Creature = std::vector<double>;

	std::vector<Creature> creatures;
	const double lower_bound;
	const double upper_bound;
	const double survival_rate;
	const double mutation_rate;
	const double crossover_rate;
	Creature bestCreature;
	double bestFitness;

	GeneticAlgorithm(const std::vector<Creature>& creatures, const double lower_bound, const double upper_bound,
					 const double mutation_rate, const double crossover_rate, const double survival_rate);
};

#endif	// GENETIC_ALGORITHM_HPP