#ifndef SIMULATED_ANNEALING_HPP
#define SIMULATED_ANNEALING_HPP

#include <vector>
#include <random> // serve per std::mt19937
#include "State.hpp"
#include "ObjectiveFunction.hpp"

class SimulatedAnnealing {
private:
    ObjectiveFunction& objective;

    int dimension;
    int maxIterations;
    int dwellIterations; // Number of iterations spent at each temperature before cooling (for thermal equilibrium)

    double initialTemperature;
    double currentTemperature;
    double temperatureScale;

    double initialStepSize;      // Initial maximum perturbation size for state changes
    double currentStepSize;      // Current maximum perturbation size (adapts during run)
    double stepSizeScale;        // Step size reduction rate (0 < scale < 1)

    double boltzmannConstant;    // Physical constant for acceptance probability

    const double lowerBound;     // Minimum allowed value in search space
    const double upperBound;     // Maximum allowed value in search space

    State bestState;             // Best solution found so far
    State currentState;          // Current working solution
    State newState;              // Candidate solution being evaluated
    
// Random engine
std::mt19937 gen;           

// Internal methods
void proposeNewSolution();   // Generates a neighboring candidate solution
void updateStepSize();       // Adapts the search step size according to scaling factor
void updateTemperature();    // Cools down the system temperature according to schedule
int equilibrate(double temperature, int iterations);  // Runs Metropolis criterion for given iterations at fixed temperature

public:
    // Costruttore unificato
    SimulatedAnnealing(ObjectiveFunction& objFunc, 
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

    double melt();                      // Heats system to find initial thermal equilibrium
    double anneal();                    // Executes cooling schedule to find optimal solution
    const State& getBestState() const;  // Returns best solution found during annealing

};

#endif // SIMULATED_ANNEALING_HPP
