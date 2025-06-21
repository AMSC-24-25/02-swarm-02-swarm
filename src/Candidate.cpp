#include <vector>
#include <limits>
#include <random>
#include <algorithm>
#include <cassert>

#include "Candidate.hpp"

// Constructor: Initializes a candidate with random values within given bounds
// and computes its initial objective value.
Candidate::Candidate(const size_t dimensions_, const double lower_bound_, const double upper_bound_, const size_t seed_,  const ObjectiveFunction& func_)
        :f0(std::numeric_limits<double>::infinity()), func(func_) {
    assert(dimensions_ > 0);
    assert(std::isfinite(lower_bound_));
    assert(std::isfinite(upper_bound_));
    assert(lower_bound_ < upper_bound_);

    std::vector<double> candidate_(dimensions_);
    std::mt19937_64 rnd{seed_};
    std::uniform_real_distribution<double> candidate_dist{lower_bound_, upper_bound_};

    for (size_t i = 0; i < dimensions_; ++i) {
        candidate_[i] = candidate_dist(rnd);
    }
    candidate = candidate_;
    f0 = func_(candidate);
};

//-----------------------------------------------------------------
// Constructor: Creates a candidate from a pre-defined vector and evaluates it.
Candidate::Candidate(const std::vector<double>& candidate_, const ObjectiveFunction &func_)
    :f0(std::numeric_limits<double>::infinity()),candidate(candidate_),  func(func_) {
    f0 = func_(candidate);
};

//-----------------------------------------------------------------
// Updates the candidate's position with a new vector and recomputes its objective value.
// This is typically called after a successful crossover.
void Candidate::updatePosition(const std::vector<double>& position) {
    candidate = position;
    f0 = func(candidate);
}




