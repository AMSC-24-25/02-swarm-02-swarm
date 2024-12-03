#include "Particle.hpp"



Particle::Particle(int dimensions, double lowerBound, double upperBound) {
    this->dimensions = dimensions;
    this->lower= lowerBound;
    this->upper = upperBound;
    position.resize(dimensions);
    velocity.resize(dimensions);
    bestLocalPosition.resize(dimensions);
    bestFitness = numeric_limits<double>::infinity();

    for (int i = 0; i < dimensions; ++i) {
        position[i] = lowerBound + static_cast<double>(rand()) / RAND_MAX * (upperBound - lowerBound);

        velocity[i] = (static_cast<double>(rand()) / RAND_MAX * 2 - 1) * (upperBound - lowerBound) * 0.1;
        bestLocalPosition[i] = position[i];
    }
}
    
void Particle::update(vector<double>& globalBestPosition) {
    for (int i = 0; i < dimensions; ++i) {
        velocity[i] =  velocity[i]*(0.4) + (0.6)*(bestLocalPosition[i]- position[i]) + (0.5) * (globalBestPosition[i] - position[i]);
        //sistemrae quando esce da lower/upper bound
        position[i] = position[i] + velocity[i];

    }
    double newValuoOfFunction = ObjectiveFunction().getValueFunction(position);
    if(newValuoOfFunction < bestFitness) {
        bestFitness = newValuoOfFunction;

        for(int i = 0; i < dimensions; ++i) {
            bestLocalPosition[i] = position[i];
        }
    }
}


