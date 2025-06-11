#ifndef CREATURE_HPP
#define CREATURE_HPP

#include <vector>

class Creature {
   public:
	std::vector<double> position;
	double fitness;

	Creature(const std::vector<double>& position);
};

#endif	// CREATURE_HPP
