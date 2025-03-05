#ifndef POSITION_HPP
#define POSITION_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>

#include "ObjectiveFunction.hpp"

class Position {
    public:
        std::vector<double> position;
        double bestFitness;
        const size_t dimensions;
        double beta;

        Position(const size_t dimensions, const double lower_bound, const double upper_bound, const size_t seed, double beta, const ObjectiveFunction& func);

        void update(const ObjectiveFunction& func, const double lower_bound,
				const double upper_bound, const size_t seed, const double gamma,
                const double sigma, const double f_thresh, const double beta_adjsut_factor);
    
};




#endif