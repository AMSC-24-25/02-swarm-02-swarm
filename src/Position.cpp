#include <vector>
#include <limits>
#include <random>
#include <algorithm>
#include <cassert>

#include "Position.hpp"

Position::Position(const size_t dimensions_, const double lower_bound, const double upper_bound, const size_t seed, const double beta_, const ObjectiveFunction& func)
    : bestFitness(std::numeric_limits<double>::infinity()), dimensions(dimensions_) {
    assert(dimensions_ > 0);
	assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
	assert(lower_bound < upper_bound);

    std::vector<double> position_(dimensions);

    std::mt19937 rnd{seed};
	std::uniform_real_distribution<double> position_dist{lower_bound, upper_bound};


    for (size_t i = 0; i < dimensions; i++) {
		position_[i] = position_dist(rnd);
	}

    position = position_;
    beta = beta_;
    bestFitness = func(position);
}


void Position::update(const ObjectiveFunction& func, const double lower_bound,
				const double upper_bound, const size_t seed, const double gamma,
                const double sigma, const double f_thresh, const double beta_adjsut_factor){


    
    
}