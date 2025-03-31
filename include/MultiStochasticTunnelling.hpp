#ifndef MULTI_STOCHASTIC_TUNNELLING_HPP
#define MULTI_STOCHASTIC_TUNNELLING_HPP

#include <vector>

#include "Position.hpp"
#include "ObjectiveFunction.hpp"

class MultiStochasticTunnelling{
    private:
    const double lower_bound;
    const double upper_bound;
    const ObjectiveFunction& func;
    const double gamma;
    const double sigma_max;
    const double sigma_min;
    const double beta_adjust_factor;
    const double max_iter;
    std::vector<std::vector<double>> candidate_positions;
    double delta;
    double beta_thresholding;


    public:
    std::vector<Position> pos;
    

    MultiStochasticTunnelling(std::vector<Position>& pos, const double lower_bound, const double upper_bound, const double sigma_max, const double sigma_min,
                        const double gamma, const double beta_adjust_factor,
                        const size_t max_iter, const ObjectiveFunction& func, double beta_thresholding);

    void iteration(const size_t seed, const size_t k);

    double mapped_function_value(const std::vector<double>&  posi);

    bool delta_condition(double delt);

    bool metropolis_condition(const double delta_f_stun, const size_t seed, const double beta, const double delta_f, const double old_delta);

    void first_k_iteration(const size_t seed, const size_t k);

    double compute_sigma(size_t i);

    void update_beta_thresholding(size_t k);
};



#endif