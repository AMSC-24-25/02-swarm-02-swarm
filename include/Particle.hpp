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
	double bestFitness;	 // Miglior fitness personale (f(bestLocalpostion)
	int dimensions;

	Particle(int dimensions, std::vector<double>& lower,
			 std::vector<double>& upper);
	void update(std::vector<double>& globalBestPosition,
				const std::vector<double>& lower,
				const std::vector<double>& upper);
};

#endif
