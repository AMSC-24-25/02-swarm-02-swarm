#ifndef SWARM_HPP
#define SWARM_HPP

#include "Particle.hpp"

class Swarm{

  public:
    vector<Particle> particles;
    vector<double> bestGlobalPosition;
    vector<double> lower;
    vector<double> upper;
    double minimus;
    Swarm(vector<Particle>& particles, vector<double> lower_, vector<double> upper_);

    double findbestFitness();

    void updateParticles();
};

#endif