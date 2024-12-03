#include <vector>

#include "Swarm.hpp"

Swarm::Swarm(std::vector<Particle>& particles_, std::vector<double> lower_,
			 std::vector<double> upper_) {
	particles = particles_;
	minimus = particles[0].bestFitness;
	bestGlobalPosition = particles[0].bestLocalPosition;
	lower = lower_;
	upper = upper_;
}

double Swarm::findbestFitness() {
	for (size_t i = 0; i < particles.size(); ++i) {
		if (particles[i].bestFitness < minimus) {
			minimus = particles[i].bestFitness;
			bestGlobalPosition = particles[i].position;
		}
	}
	return minimus;
}

void Swarm::updateParticles() {
	for (size_t i = 0; i < particles.size(); ++i) {
		particles[i].update(bestGlobalPosition);
	}
}
