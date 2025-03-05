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

#include "StochasticTunnelling.hpp"

StochasticTunnelling::StochasticTunnelling(const Position& pos_, const double lower_bound_, const double upper_bound_,
                        const double gamma_, const double sigma_, const double f_thresh_, const double beta_adjust_factor_,
                        const double max_iter_, ObjectiveFunction& func_)
    : lower_bound(lower_bound_),
      upper_bound(upper_bound_),
      func(func_),
      pos(pos),
      gamma(gamma_),
      sigma(sigma_),
      f_thresh(f_thresh_),
      beta_adjust_factor(beta_adjust_factor_),
      max_iter(max_iter_) {
    assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
    assert(lower_bound < upper_bound);
    assert(beta_adjust_factor > 0.0 && beta_adjust_factor < 1.0);
    assert(f_thresh > 0.0 && f_thresh < 1.0);
}


void StochasticTunnelling::iteration(){}

