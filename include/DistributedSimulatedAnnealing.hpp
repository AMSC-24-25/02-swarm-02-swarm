#ifndef DISTRIBUTED_SIMULATED_ANNEALING_HPP
#define DISTRIBUTED_SIMULATED_ANNEALING_HPP

#include <vector>
#include <random> // serve per std::mt19937
#include "State.hpp"
#include "ObjectiveFunction.hpp"

class DistributedSimulatedAnnealing {
private:
    ObjectiveFunction& objective;

    int dimension;
    int maxIterations;
    int dwellIterations;

    double initialTemperature;
    double currentTemperature;
    double temperatureScale;

    double initialStepSize;
    double currentStepSize;
    double stepSizeScale;

    double boltzmannConstant;

    

    const double lowerBound;
    const double upperBound;


   
   
    State bestState;
    State currentState;
    State newState;
    
    // Random engine
    std::mt19937 gen;

    // Metodi interni
    void proposeNewSolution();
    void updateStepSize();
    void updateTemperature();
    int equilibrate(double temperature, int iterations);

public:
    // Costruttore unificato
    DistributedSimulatedAnnealing(ObjectiveFunction& objFunc, 
                       int dim,
                       int maxIter, 
                       int dwell,
                       double temp, 
                       double tempScale,
                       double stepSize, 
                       double stepScale,
                       double boltzmannK,
                       const std::vector<double>& initialGuess,
                       const double lowerBound,
                       const double upperBound,
                       size_t seed 
                      );

    double melt();
    double anneal();
    const State& getBestState() const; 

};

#endif // SIMULATED_ANNEALING_HPP
