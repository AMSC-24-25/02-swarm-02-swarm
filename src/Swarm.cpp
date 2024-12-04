#include <vector>

#include "Swarm.hpp"

Swarm::Swarm(std::vector<Particle>& particles_, std::vector<double> lower_, std::vector<double> upper_ , double c1_, double c2_ ) { 
    c1=c1_;
	c2=c2_;
	particles = particles_;
	minimum = particles[0].bestFitness;
	bestGlobalPosition = particles[0].bestLocalPosition;
	lower = lower_;
	upper = upper_;
}

double Swarm::findbestFitness() {
	for (size_t i = 0; i < particles.size(); ++i) {
		if (particles[i].bestFitness < minimum) {
			minimum = particles[i].bestFitness;
			bestGlobalPosition = particles[i].position;
		}
	}
	return minimum;
}

void Swarm::updateParticles(double c1, double c2, double w) {
	for (size_t i = 0; i < particles.size(); ++i) {
		particles[i].update(bestGlobalPosition, lower, upper, c1, c2, w);
	}
}
