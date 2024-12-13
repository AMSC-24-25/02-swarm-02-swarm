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
	double minimum;
	double c1;
	double c2;
	double w;

	Swarm(std::vector<Particle>& particles, std::vector<double> lower_, std::vector<double> upper_, double c1_,
		  double c2_, double w_);

	double findbestFitness();

	void updateParticles();

	void updateInertia(int max_iterations, double w_min, double w_max);
};

#endif
