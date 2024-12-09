#include <vector>

#include "Swarm.hpp"

Swarm::Swarm(const std::vector<Particle>& particles_, const std::vector<double> lower_,
			 const std::vector<double> upper_, const double c1_, const double c2_, const double w_,
			 ObjectiveFunction& func_)
	: particles(particles_),
	  bestGlobalPosition(particles[0].bestLocalPosition),
	  lower(lower_),
	  upper(upper_),
	  minimum(particles[0].bestFitness),
	  c1(c1_),
	  c2(c2_),
	  w(w_),
	  func(func_) {}

double Swarm::findBestFitness() {
	for (size_t i = 0; i < particles.size(); ++i) {
		if (particles[i].bestFitness < minimum) {
			minimum = particles[i].bestFitness;
			bestGlobalPosition = particles[i].position;
		}
	}
	return minimum;
}

void Swarm::updateParticles() {
	for (size_t i = 0; i < particles.size(); ++i) {
		particles[i].update(func, bestGlobalPosition, lower, upper, c1, c2, w);
	}
}

void Swarm::updateInertia(const int max_iterations, const double w_min, const double w_max) {
	w = w - ((w_max - w_min) / max_iterations);
}
