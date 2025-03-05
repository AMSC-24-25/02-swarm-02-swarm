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


void StochasticTunnelling::iteration(const size_t seed){
  candidate_position = pos.generate_new_position(func, lower_bound, upper_bound, seed, sigma);

  delta = mapped_function_value(candidate_position) - mapped_function_value(pos.position);

  if(delta_condition(mapped_function_value(candidate_position), mapped_function_value(pos.position)) || metropolis_condition(delta, seed, pos.beta)){
    pos.update_position(candidate_position, func);
    
    pos.update_avg_window(mapped_function_value(pos.position));

    pos.update_beta(f_thresh, beta_adjust_factor);
  }
  

  
}


double StochasticTunnelling::mapped_function_value(const std::vector<double>& posi){
  
  return (1.0 + std::exp(-gamma * (func(posi) - pos.f0)));

}

bool StochasticTunnelling::delta_condition(double map_value_new, double map_value_old){
  double delta = map_value_new - map_value_old;
  if(delta <= 0){
    return true;
  }else{
    return false;
  }
}

bool StochasticTunnelling::metropolis_condition(const double delta_f_stun, const size_t seed, const double beta) {
        std::mt19937 gen(seed); 
        
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        double random_value = dist(gen);

        double exp_value = std::exp(-beta * delta_f_stun);

        return random_value < exp_value;
}

