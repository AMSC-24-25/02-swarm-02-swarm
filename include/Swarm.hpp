#ifndef SWARM_HPP
#define SWARM_HPP

#include <vector>

#include "Particle.hpp"

class Swarm {
   public:
	std::vector<Particle> particles;
	std::vector<double> bestGlobalPosition;
	std::vector<double> lower;
	std::vector<double> upper;
	double minimus;
	Swarm(std::vector<Particle>& particles, std::vector<double> lower_,
		  std::vector<double> upper_);

	double findbestFitness();

	void updateParticles();
};

#endif
