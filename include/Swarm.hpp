#ifndef SWARM_HPP
#define SWARM_HPP

#include <vector>

#include "Particle.hpp"
#include "ObjectiveFunction.hpp"

class Swarm {
   public:
	std::vector<Particle> particles;
	std::vector<double> bestGlobalPosition;
	double minimum;
	const double lower_bound;
	const double upper_bound;
	const double c1;
	const double c2;
	double w;
	size_t seed;
	const size_t n_threads;
	const ObjectiveFunction& func;

	Swarm(const std::vector<Particle>& particles, const double lower_bound_, const double upper_bound_,
		  const double c1_, const double c2_, const double w_, const size_t seed, ObjectiveFunction& func_,
		  const size_t n_threads);

	void findBestFitness();

	void updateParticles();

	void updateInertia(const size_t max_iterations, const double w_min, const double w_max);
};

#endif
