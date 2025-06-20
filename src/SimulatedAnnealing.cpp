#include "SimulatedAnnealing.hpp"
#include "State.hpp"
#include <random>
#include <cmath>
#include <limits>
#include <iostream>
#include <cassert>
#include <algorithm>

SimulatedAnnealing::SimulatedAnnealing(ObjectiveFunction& objFunc,
                                       int dim,
                                       int maxIter,
                                       int dwell, 
                                       double temp,
                                       double tempScale,
                                       double stepSize,
                                       double stepScale,
                                       double boltzmannK,
                                       const std::vector<double>& initialGuess,
                                       double lowerBound_,
                                       double upperBound_,
                                       size_t seed)
    : objective(objFunc),
      dimension(dim),
      maxIterations(maxIter),
      dwellIterations(dwell),
      initialTemperature(temp),
      currentTemperature(temp),
      temperatureScale(tempScale),
      initialStepSize(stepSize),
      currentStepSize(stepSize),
      stepSizeScale(stepScale),
      boltzmannConstant(boltzmannK),
      lowerBound(lowerBound_),
      upperBound(upperBound_),
      bestState(initialGuess, objFunc),
      currentState(initialGuess, objFunc),
      newState(std::vector<double>(dim, 0.0), objFunc),
      gen(seed)

     
{   assert(std::isfinite(lowerBound));
    assert(std::isfinite(upperBound));
    assert(lowerBound < upperBound);

    for (int i = 0; i < dimension; ++i) {
        assert(initialGuess[i] >= lowerBound && initialGuess[i] <= upperBound);
    }
}

// Generates a new candidate solution by randomly perturbing the current state
// Adds uniformly distributed random values (range: ±currentStepSize) to each dimension
// Maintains search within specified bounds through clamping

void SimulatedAnnealing::proposeNewSolution() {
    std::uniform_real_distribution<> dist(-currentStepSize, currentStepSize);
    std::vector<double> candidateValues(dimension);
    for (int i = 0; i < dimension; ++i) {
        double candidate = currentState.values[i] + dist(gen);
        candidateValues[i] = std::clamp(candidate, lowerBound, upperBound);
    }
    newState.update(candidateValues);
}

// Metropolis acceptance: 
// - Always accepts improving moves (ΔE < 0)
// - Accepts degrading moves with P = exp(-ΔE/(k·T)) where:
//   ΔE = energy difference, k = Boltzmann constant, T = current temperature
// compare against uniform random [0,1] for probabilistic acceptance
        
int SimulatedAnnealing::equilibrate(double temperature, int iterations) {
    int accepted = 0;
    for (int i = 0; i < iterations; ++i) {
        proposeNewSolution();
        double newCost = newState.cost;
        double delta = newCost - currentState.cost;

        if (delta < 0 || std::exp(-delta / (boltzmannConstant * temperature)) > std::generate_canonical<double, 10>(gen)) { //Genera un numero casuale tra 0 e 1 con alta precisione (10 bit di randomicità).
            
            currentState.update(newState.values);

            if (newCost < bestState.cost) {
                bestState.update(newState.values);

            }
            ++accepted;
        }
    }
    return accepted;
}


// Heating phase - raises temperature until:
// 1. Acceptance ratio drops below 10% of dwell iterations
// 2. Temperature increases linearly (+1.0 per iteration)
// Updates and returns final pre-annealing temperature


double SimulatedAnnealing::melt() {
    int i = 0;
    double temp = currentTemperature;
    while (i < maxIterations) {
        int accepted = equilibrate(temp, dwellIterations);
        if (accepted < dwellIterations / 10) break;
        temp += 1.0;
        ++i;
    }
    currentTemperature = temp;
    return temp;
}

// Annealing process: gradually cools system by multiplying temperature by scaling factor

double SimulatedAnnealing::anneal() {
for (int i = 0; i < maxIterations; ++i) {
    equilibrate(currentTemperature, dwellIterations);
    updateTemperature();
    updateStepSize();
    

    // Se la temperatura corrente è molto bassa, interrompiamo l'annealing
    if (currentTemperature < 1e-6) break;
}
return currentTemperature;
}


// Updates current temperature using the temperature scaling factor (geometric cooling)

void SimulatedAnnealing::updateTemperature() {
    currentTemperature *= temperatureScale;
}

// Updates step size using the scaling factor (geometric reduction for finer search)

void SimulatedAnnealing::updateStepSize() {
    currentStepSize *= stepSizeScale;
}

const State& SimulatedAnnealing::getBestState() const {
    return bestState;
}


