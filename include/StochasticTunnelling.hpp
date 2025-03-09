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
    const double sigma_max;
    const double sigma_min;
    const double f_thresh;
    const double beta_adjust_factor;
    const double max_iter;
    std::vector<double> candidate_position;
    double delta;


    public:
    Position pos;
    

    StochasticTunnelling(Position& pos, const double lower_bound, const double upper_bound, const double sigma_max, const double sigma_min,
                        const double gamma, const double f_thresh, const double beta_adjust_factor,
                        const size_t max_iter, const ObjectiveFunction& func);

    void iteration(const size_t seed, const size_t k);

    double mapped_function_value(const std::vector<double>&  posi);

    bool delta_condition(double delt);

    bool metropolis_condition(const double delta_f_stun, const size_t seed, const double beta);

    void first_k_iteration(const size_t seed, const size_t k);

    double compute_sigma(size_t i);
};



#endif // STOCHASTIC_TUNNELLING_HPP