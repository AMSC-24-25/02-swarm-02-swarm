#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>

#include "ObjectiveFunction.hpp"

class Particle {
   public:
	std::vector<double> position;			// Posizione attuale
	std::vector<double> velocity;			// Velocità attuale
	std::vector<double> bestLocalPosition;	// Miglior posizione personale
	double bestFitness;						// Miglior fitness personale (f(bestLocalPosition)
	int dimensions;

	Particle(const int dimensions, const std::vector<double>& lower, const std::vector<double>& upper);

	void update(const ObjectiveFunction& func, const std::vector<double>& globalBestPosition,
				const std::vector<double>& lower, const std::vector<double>& upper, const double c1, const double c2,
				const double w);
};

#endif
