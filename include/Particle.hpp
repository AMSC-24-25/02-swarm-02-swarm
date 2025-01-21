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
	const size_t dimensions;

	Particle(const size_t dimensions, const double lower_bound, const double upper_bound, const size_t seed);

	void update(const ObjectiveFunction& func, const std::vector<double>& globalBestPosition, const double lower_bound,
				const double upper_bound, const double c1, const double c2, const double w, const size_t seed);
};

#endif
