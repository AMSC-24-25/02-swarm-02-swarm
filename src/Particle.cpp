#include "Particle.hpp"



Particle::Particle(int dimensions, vector<double> lowerBound, vector<double> upperBound) {
    this->dimensions = dimensions;
    lower= lowerBound;
    upper = upperBound;
    std::cout<<lower.size()<<" "<<upper.size()<<" "<<dimensions<<std::endl;
    bestFitness = numeric_limits<double>::infinity();

    for (int i = 0; i < dimensions; ++i) {
        position.push_back(lower[i] + static_cast<double>(rand()) / RAND_MAX * (upper[i] - lower[i]));
        velocity.push_back((static_cast<double>(rand()) / RAND_MAX * 2 - 1) * (upper[i] - lower[i]) * 0.1);
        bestLocalPosition.push_back(position[i]);
        std::cout<<position[i]<<" "<<velocity[i]<<std::endl;
        std::cout<<upper[i]<<" "<<lower[i]<<std::endl;
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


