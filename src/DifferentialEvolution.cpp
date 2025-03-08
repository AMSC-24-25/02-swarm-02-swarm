#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>
#include <limits>
#include <utility>
#include <numeric>
#include <omp.h>

#include "DifferentialEvolution.hpp"

#include <memory>

DifferentialEvolution::DifferentialEvolution(const std::vector<Candidate>& candidates_, const size_t dimensions_, const double lower_bound_, const double upper_bound_, const size_t seed_, const size_t max_gen_, const double F_, const double CR_, ObjectiveFunction& func_, const size_t n_threads_)
    : candidates(candidates_), dimensions(dimensions_),lower_bound(lower_bound_),upper_bound(upper_bound_), seed(seed_),F(F_), CR(CR_),max_gen(max_gen_) , func(func_), n_threads(n_threads_) {};

void DifferentialEvolution::select_three_random(int excluded_index,int& i1, int& i2, int& i3, std::mt19937 &gen) {
    std::uniform_int_distribution<int> dist(0,candidates.size()-1);
    do {
        i1 = dist(gen);
    }while (i1 == excluded_index);

    do {
        i2 = dist(gen);
    }while (i2 == excluded_index || i2 == i1);

    do {
        i3 = dist(gen);
    }while (i3 == excluded_index || i3 == i1 || i3 == i2);
};

void DifferentialEvolution::mutate(std::vector<double>& mutantPos, const Candidate &r1, const Candidate &r2, const Candidate &r3) {
    for (int i = 0; i < dimensions; i++) {
        mutantPos[i] = r1.candidate[i] + F * (r2.candidate[i] - r3.candidate[i]);
    }
};

void DifferentialEvolution::crossover(Candidate& original, Candidate mutant, std::mt19937 &gen) {
    std::vector<double> mergedPos= original.candidate;
    std::uniform_int_distribution<int> intDist(0,dimensions-1);
    int j_rand = intDist(gen);

    std::uniform_real_distribution<double> realDist(0.0,1.0);
    for (size_t j = 0; j< dimensions; j++) {
        if (realDist(gen)<CR || j == j_rand) {
            mergedPos[j] = mutant.candidate[j];
        }
    }
    std::shared_ptr<Candidate> merged= std::make_shared<Candidate>(mergedPos,func);
    if (merged->f0 < original.f0) {
        original.updatePosition(mergedPos);
    }
}

void DifferentialEvolution::updateCandidate() {
    std::random_device rd;
    std::mt19937 gen(rd());

    for (size_t i = 0; i < candidates.size(); i++) {
        int i1;
        int i2;
        int i3;
        bool choice;
        select_three_random(i,i1,i2,i3,gen);
        Candidate r1 = candidates[i1];
        Candidate r2 = candidates[i2];
        Candidate r3 = candidates[i3];

        std::vector<double> mutantPos;
        mutantPos.resize(dimensions);
        mutate(mutantPos,r1,r2,r3);
        std::shared_ptr<Candidate> mutant= std::make_shared<Candidate>(mutantPos,func);

        crossover(candidates[i], (*mutant),gen);
    }
};

int DifferentialEvolution::findSol() {
    for (int i = 0; i < candidates.size(); i++) {
        if (candidates[i].f0 < 1e-5) {
            return i;
        }
    }
    return -1;
}

void DifferentialEvolution::iterate() {
    int solPos;
    int iteration=0;
    while (iteration<max_gen) {
        iteration++;
        solPos = findSol();
        if (solPos != -1) {
            break;
        }
        updateCandidate();
    }
}

