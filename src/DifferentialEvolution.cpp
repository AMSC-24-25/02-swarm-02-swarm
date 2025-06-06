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

// Constructor for DifferentialEvolution.
// Initializes key parameters and verifies that the bounds and control parameters are valid.
DifferentialEvolution::DifferentialEvolution(const std::vector<Candidate>& candidates_, const size_t dimensions_, const double lower_bound_, const double upper_bound_, const size_t seed_, const size_t max_gen_, const double F_, const double CR_, ObjectiveFunction& func_, const size_t n_threads_)
    : candidates(candidates_), dimensions(dimensions_),lower_bound(lower_bound_),upper_bound(upper_bound_), seed(seed_),F(F_), CR(CR_),max_gen(max_gen_) , func(func_), n_threads(n_threads_),gens(n_threads_), bestCandidate(nullptr) {
    assert(candidates.size() >= 4 && "Differential Evolution requires at least 4 candidates");
    assert(std::isfinite(lower_bound));
    assert(std::isfinite(upper_bound));
    assert(lower_bound < upper_bound);
    assert(F>=0.0 && F <=1.0);
    assert(CR>=0.0 && CR<=1.0);

    for (size_t t = 0; t < n_threads; ++t) {
        gens[t].seed(seed + t);
    }

};

//-----------------------------------------------------------------
// Selects three distinct random indices (different from the target index)
// from the candidate pool for mutation operations.
void DifferentialEvolution::select_three_random(int excluded_index,int& i1, int& i2, int& i3, std::mt19937& local_gen) {
    std::uniform_int_distribution<int> dist(0,candidates.size()-1);
    do {
        i1 = dist(local_gen);
    }while (i1 == excluded_index);

    do {
        i2 = dist(local_gen);
    }while (i2 == excluded_index || i2 == i1);

    do {
        i3 = dist(local_gen);
    }while (i3 == excluded_index || i3 == i1 || i3 == i2);
};

//-----------------------------------------------------------------
// Generates a mutant vector using the DE/rand/1 strategy.
// Each component is computed as: r1 + F * (r2 - r3) and clamped within bounds.
void DifferentialEvolution::mutate(std::vector<double>& mutantPos, const Candidate &r1, const Candidate &r2, const Candidate &r3) {
    for (size_t i = 0; i < dimensions; i++) {
        mutantPos[i] =std::clamp(r1.candidate[i] + F * (r2.candidate[i] - r3.candidate[i]), lower_bound, upper_bound);
    }
};

//-----------------------------------------------------------------
// Performs crossover between the original candidate and the mutant.
// Ensures that at least one dimension (j_rand) comes from the mutant vector.
// If the trial candidate has a better objective value, it replaces the original.
void DifferentialEvolution::crossover(Candidate& original,const Candidate& mutant, std::mt19937& local_gen) {
    std::vector<double> mergedPos= original.candidate;
    std::uniform_int_distribution<int> intDist(0,dimensions-1);
    size_t j_rand = intDist(local_gen); // Guarantee at least one mutant gene is selected

    std::uniform_real_distribution<double> realDist(0.0,1.0);
    for (size_t j = 0; j< dimensions; j++) {
        if (realDist(local_gen)<CR || j == j_rand) {
            mergedPos[j] = mutant.candidate[j];
        }
    }
    Candidate merged(mergedPos, func);
    if (merged.f0 < original.f0) {
        original.updatePosition(mergedPos);
    }
}

//-----------------------------------------------------------------
// Evolves each candidate in the population in parallel using OpenMP.
// For every candidate, three random distinct candidates are chosen for mutation,
// followed by the mutation and crossover steps to generate a trial candidate.
void DifferentialEvolution::updateCandidate() {
#pragma omp parallel num_threads(n_threads)
    {
        int tid = omp_get_thread_num();
        auto& local_gen = gens[tid];
#pragma omp for schedule(static)
        for (size_t i = 0; i < candidates.size(); i++) {
            int i1;
            int i2;
            int i3;

            select_three_random(i, i1, i2, i3, local_gen);
            Candidate r1 = candidates[i1];
            Candidate r2 = candidates[i2];
            Candidate r3 = candidates[i3];

            std::vector<double> mutantPos;
            mutantPos.resize(dimensions);
            mutate(mutantPos, r1, r2, r3);
            Candidate mutant(mutantPos, func);

            crossover(candidates[i], mutant, local_gen);
        }
    }

};

//-----------------------------------------------------------------
// Searches for the candidate with the best (minimum) objective value
// and stores a reference to it
void DifferentialEvolution::findSol() {

    double min = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < candidates.size(); i++) {

        if (candidates[i].f0 < min) {
            bestCandidate = &candidates[i];
            min = bestCandidate->f0;
        }
    }
}

