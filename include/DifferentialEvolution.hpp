#ifndef DIFFERENTIAL_EVOLUTION_HPP
#define DIFFERENTIAL_EVOLUTION_HPP

#include <vector>

#include "Candidate.hpp"
#include "ObjectiveFunction.hpp"

class DifferentialEvolution{

public:

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
    std::vector<std::mt19937> gens;


    void select_three_random(int excluded_index,int& i1, int& i2, int& i3, std::mt19937& local_gen) const;
    void mutate(std::vector<double>& mutantPos, const Candidate& r1, const Candidate& r2, const Candidate& r3) const;
    void crossover(Candidate& original,const Candidate& mutant,std::mt19937& local_gen) const;


    DifferentialEvolution(const std::vector<Candidate>& candidates_,size_t dimensions_, double lower_bound_, double upper_bound_, size_t seed_, size_t max_gen_, double F_, double CR_, ObjectiveFunction& func_, size_t n_threads_) noexcept;
    void updateCandidate();
    void findSol() noexcept;
    const Candidate& getBestCandidate() const noexcept;
    std::size_t bestIndex{0};

};



#endif // DIFFERENTIAL_EVOLUTION_HPP