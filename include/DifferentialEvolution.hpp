#ifndef DIFFERENTIAL_EVOLUTION_HPP
#define DIFFERENTIAL_EVOLUTION_HPP

#include <vector>

#include "Candidate.hpp"
#include "ObjectiveFunction.hpp"

class DifferentialEvolution{

    private:
    std::vector<Candidate> candidates;
    const size_t dimensions;
    const double lower_bound;
    const double upper_bound;
    const size_t seed;
    const double F;
    const double CR;
    const size_t max_gen;
    const ObjectiveFunction& func;
    size_t n_threads;
    std::mt19937 gen;


    void select_three_random(int excluded_index,int& i1, int& i2, int& i3);
    void mutate(std::vector<double>& mutantPos, const Candidate& r1, const Candidate& r2, const Candidate& r3);
    void crossover(Candidate& original, Candidate& mutant);


    public:
    DifferentialEvolution(const std::vector<Candidate>& candidates,const size_t dimensions_, const double lower_bound_, const double upper_bound_, const size_t seed_, const size_t max_gen_, const double F_, const double CR_, ObjectiveFunction& func_, const size_t n_threads_);
    void updateCandidate();
    void findSol();
    Candidate* bestCandidate;
};



#endif // DIFFERENTIAL_EVOLUTION_HPP