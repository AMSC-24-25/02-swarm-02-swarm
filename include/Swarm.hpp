#ifndef SWARM_HPP
#define SWARM_HPP

#include <vector>

#include "Particle.hpp"
#include "ObjectiveFunction.hpp"

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
	const size_t seed;
	const size_t n_threads;
	ObjectiveFunction& func;

	Swarm(const std::vector<Particle>& particles, const std::vector<double>& lower_, const std::vector<double>& upper_,
		  const double c1_, const double c2_, const double w_, const size_t seed, ObjectiveFunction& func_,
		  const size_t n_threads);

	double findBestFitness();

	void updateParticles();

	void updateInertia(const int max_iterations, const double w_min, const double w_max);
};

#endif
