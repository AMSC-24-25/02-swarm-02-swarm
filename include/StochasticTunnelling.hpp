#ifndef STOCHASTIC_TUNNELLING_HPP
#define STOCHASTIC_TUNNELLING_HPP

#include <vector>

#include "Position.hpp"
#include "ObjectiveFunction.hpp"

class StochasticTunnelling{
    private:
    const double lower_bound;
    const double upper_bound;
    const ObjectiveFunction& func;
    const double gamma;
    const double sigma;
    const double f_thresh;
    const double beta_adjust_factor;
    const double max_iter;
    std::vector<double> candidate_position;
    double delta;


    public:
    Position pos;
    

    StochasticTunnelling(const Position& pos, const double lower_bound, const double upper_bound,
                        const double gamma, const double sigma, const double f_thresh, const double beta_adjust_factor,
                        const double max_iter, ObjectiveFunction& func);

    void iteration(const size_t seed);

    double mapped_function_value(const std::vector<double>&  posi);

    bool delta_condition(double map_value_new, double map_value_old);

    bool metropolis_condition(const double delta_f_stun, const size_t seed, const double beta);
};



#endif // STOCHASTIC_TUNNELLING_HPP