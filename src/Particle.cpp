#include "Particle.hpp"



Particle::Particle(int dimensions_, vector<double>& lower, vector<double>& upper) {
    dimensions = dimensions_;
    bestFitness = numeric_limits<double>::infinity();
    vector<double> position_(dimensions);
    vector<double> velocity_(dimensions);
    vector<double> bestLocalPosition_(dimensions);


    for (int i = 0; i < dimensions; ++i) {
        position_[i] = (lower[i] + static_cast<double>(rand()) / RAND_MAX * (upper[i] - lower[i]));
        velocity_[i] = ((static_cast<double>(rand()) / RAND_MAX * 2 - 1) * (upper[i] - lower[i]) * 0.1);
        bestLocalPosition_[i] = (position_[i]);
    }
    position = position_;
    velocity = velocity_;
    bestLocalPosition = bestLocalPosition_;
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


