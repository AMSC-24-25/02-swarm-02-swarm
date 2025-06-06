#include "SimulatedAnnealing.hpp"
#include "State.hpp"
#include <random>
#include <cmath>
#include <limits>
#include <iostream>
#include <cassert>
#include <algorithm>


//Costruttore: inizializza i parametri dell'algoritmo

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

     
{
    assert(initialGuess.size() == dimension);
    assert(std::isfinite(lowerBound));
    assert(std::isfinite(upperBound));
    assert(lowerBound < upperBound);

    for (int i = 0; i < dimension; ++i) {
        assert(initialGuess[i] >= lowerBound && initialGuess[i] <= upperBound);
    }
}



// Proposta di una nuova soluzione: genera una nuova soluzione perturbata in base alla soluzione corrente e alla dimensione del passo
// La nuova soluzione è generata aggiungendo un valore casuale alla soluzione corrente
// La distribuzione uniforme è utilizzata per generare un numero casuale compreso tra -currentStepSize e currentStepSize

void SimulatedAnnealing::proposeNewSolution() {
    std::uniform_real_distribution<> dist(-currentStepSize, currentStepSize);
    std::vector<double> candidateValues(dimension);
    for (int i = 0; i < dimension; ++i) {
        double candidate = currentState.values[i] + dist(gen);
        candidateValues[i] = std::clamp(candidate, lowerBound, upperBound);
    }
    newState.update(candidateValues);
}

// Equilibrio: accetta o rifiuta la nuova soluzione in base alla probabilità di accettazione (Metropolis criterion). Se la nuova soluzione è migliore, viene accettata. 
//Se è peggiore, viene accettata con una certa probabilità in base alla temperatura corrente e alla costante di Boltzmann
//Iterare un po' di volte per ottenere l'equilibrio termico
// Se la nuova soluzione è migliore, accettala, altrimenti accettala con una certa probabilità pari a exp(-delta / (k * temperature)) 
// std::generate_canonical<double, 10>(gen) genera un numero casuale tra 0 e 1 e serve per confrontare con la probabilità. Se la probabilità è più alta 
//del numero casuale → accettiamo comunque la soluzione peggiore perchè all'inizio (temperatura alta), l’algoritmo esplora più liberamente, mentre alla fine 
//(temperatura bassa), si concentra sulle soluzioni migliori così si evitano minimi locali.
        
        
int SimulatedAnnealing::equilibrate(double temperature, int iterations) {
    int accepted = 0;
    for (int i = 0; i < iterations; ++i) {
        proposeNewSolution();
        double newCost = newState.cost;
        double delta = newCost - currentState.cost;

        if (delta < 0 || std::exp(-delta / (boltzmannConstant * temperature)) > std::generate_canonical<double, 10>(gen)) {
            
            currentState.update(newState.values);

            if (newCost < bestState.cost) {
                bestState.update(newState.values);

            }
            ++accepted;
        }
    }
    return accepted;
}


// Funzione di fusione: aumenta la temperatura corrente fino a quando il sistema "si scioglie" (il numero di soluzioni accettate è sufficientemente alto)
// La temperatura viene aumentata di 1.0 ad ogni iterazione fino a quando il numero di soluzioni accettate è inferiore a un decimo delle iterazioni di equilibrio
// La temperatura corrente viene aggiornata alla fine della fusione
// La funzione restituisce la temperatura finale
/* equilibrium is defined as a 10% or smaller change in 10 iterations */


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

//Raffreddamento: riduce la temperatura corrente in base al fattore di scala della temperatura

double SimulatedAnnealing::anneal() {
    for (int i = 0; i < maxIterations; ++i) {
        equilibrate(currentTemperature, dwellIterations);
        updateTemperature();
        updateStepSize();


    }
    return currentTemperature;
}

// Funzione per aggiornare la temperatura corrente in base al fattore di scala della temperatura

void SimulatedAnnealing::updateTemperature() {
    currentTemperature *= temperatureScale;
}

// Funzione per aggiornare la dimensione del passo corrente in base al fattore di scala della dimensione del passo
// La dimensione del passo viene ridotta ad ogni iterazione per consentire una ricerca più fine delle soluzioni migliori


void SimulatedAnnealing::updateStepSize() {
    currentStepSize *= stepSizeScale;
}

const State& SimulatedAnnealing::getBestState() const {
    return bestState;
}


