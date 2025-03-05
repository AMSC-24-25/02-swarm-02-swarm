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


    public:
    Position pos;
    

    StochasticTunnelling(const Position& pos, const double lower_bound, const double upper_bound,
                        const double gamma, const double sigma, const double f_thresh, const double beta_adjust_factor,
                        const double max_iter, ObjectiveFunction& func);

    void iteration();

};



#endif // STOCHASTIC_TUNNELLING_HPP