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

StochasticTunnelling::StochasticTunnelling(Position& pos_, const double lower_bound_, const double upper_bound_,const double sigma_max_, const double sigma_min_,
                        const double gamma_, const double beta_adjust_factor_,
                        const size_t max_iter_, const ObjectiveFunction& func_)
    : lower_bound(lower_bound_),
      upper_bound(upper_bound_),
      func(func_),
      gamma(gamma_),
      sigma_max(sigma_max_),
      sigma_min(sigma_min_),
      beta_adjust_factor(beta_adjust_factor_),
      max_iter(max_iter_),
      pos(pos_) {
    assert(std::isfinite(lower_bound));
	  assert(std::isfinite(upper_bound));
    assert(lower_bound < upper_bound);
    assert(beta_adjust_factor > 0.0 && beta_adjust_factor < 1.0);
}


void StochasticTunnelling::iteration(const size_t seed, const size_t k){
  double sigma = compute_sigma(k);
  candidate_position = pos.generate_new_position(lower_bound, upper_bound, seed, sigma);
  std::cout<<"current position: "<<pos.position[0]<<" "<<pos.position[1]<<std::endl;
  std::cout<<"value func: "<<func(pos.position)<<" value mapped: "<<mapped_function_value(pos.position)<<std::endl;
  std::cout<<"new position: "<<candidate_position[0]<<" "<<candidate_position[1]<<std::endl;
  std::cout<<"new value func: "<<func(candidate_position)<<" new value mapped: "<<mapped_function_value(candidate_position)<<std::endl;
  

  delta = mapped_function_value(candidate_position) - mapped_function_value(pos.position);
  std::cout<<"delta function mapped: "<<delta<<std::endl;


  //if(delta_condition(delta) || metropolis_condition(delta, seed, pos.beta)){
  if(delta_condition(delta) or metropolis_condition(delta, seed, pos.beta)){
    pos.update_position(candidate_position, func);
    
    pos.update_avg_window(mapped_function_value(pos.position));

  }

  pos.update_beta(beta_adjust_factor);
  
}


void StochasticTunnelling::first_k_iteration(const size_t seed, const size_t k){
  double sigma = compute_sigma(k);
  candidate_position = pos.generate_new_position(lower_bound, upper_bound, seed, sigma);
  std::cout<<"current position: "<<pos.position[0]<<" "<<pos.position[1]<<std::endl;
  std::cout<<"value func: "<<func(pos.position)<<" value mapped: "<<mapped_function_value(pos.position)<<std::endl;
  std::cout<<"new position: "<<candidate_position[0]<<" "<<candidate_position[1]<<std::endl;
  std::cout<<"new value func: "<<func(candidate_position)<<" new value mapped: "<<mapped_function_value(candidate_position)<<std::endl;


  delta = mapped_function_value(candidate_position) - mapped_function_value(pos.position);
  std::cout<<"delta function mapped: "<<delta<<std::endl;

  if(delta_condition(delta) or metropolis_condition(delta, seed, pos.beta)){
    pos.update_position(candidate_position, func);
    
    pos.increase_avg_window_at_position(mapped_function_value(pos.position), k);

  }
}


double StochasticTunnelling::mapped_function_value(const std::vector<double>& posi){  
  return (1.0 - std::exp(-gamma * (func(posi) - pos.f0)));
}

bool StochasticTunnelling::delta_condition(double delt){
  if(delt <= 0){
    return true;
  }else{
    return false;
  }
}

bool StochasticTunnelling::metropolis_condition(const double delta_f_stun, const size_t seed, const double beta) {
        static std::mt19937 gen(seed); 
        
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        double random_value = dist(gen);

        double exp_value = std::exp(-beta * delta_f_stun);

        std::cout<<"beta: "<<beta<<std::endl;

        std::cout<<"percentage: "<<exp_value<<std::endl;

        return random_value < exp_value;
}

double StochasticTunnelling::compute_sigma(size_t i){
  return sigma_max - (((sigma_max - sigma_min)/ max_iter))*i;
}
