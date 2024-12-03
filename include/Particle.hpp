#ifndef PARTICLE_HPP
#define PARTICLE_HPP


#include <bits/stdc++.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>
#include "ObjectiveFunction.hpp"


class Particle {
public:
    vector<double> position;    // Posizione attuale
    vector<double> velocity;    // Velocità attuale
    vector<double> bestLocalPosition; // Miglior posizione personale
    double bestFitness;         // Miglior fitness personale (f(bestLocalpostion)
    vector<double> lower;
    vector<double> upper;
    int dimensions;

    Particle(int dimensions, vector<double> lowerBound, vector<double> upperBound);
    void update(vector<double>& globalBestPosition);


};

#endif