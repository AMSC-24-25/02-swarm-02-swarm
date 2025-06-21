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

#include "MultiStochasticTunnelling.hpp"

MultiStochasticTunnelling::MultiStochasticTunnelling(std::vector<Position>& pos_, const double lower_bound_, const double upper_bound_,const double sigma_max_, const double sigma_min_,
                        const double gamma_, const double beta_adjust_factor_,
                        const size_t max_iter_, const ObjectiveFunction& func_, double beta_thresholding_, const size_t num_positions_, const size_t time_step_updating_, const size_t dimension_, const size_t n_threads_)
    : lower_bound(lower_bound_),
      upper_bound(upper_bound_),
      func(func_),
      gamma(gamma_),
      sigma_max(sigma_max_),
      sigma_min(sigma_min_),
      beta_adjust_factor(beta_adjust_factor_),
      max_iter(max_iter_),
      beta_thresholding(beta_thresholding_),
      num_positions(num_positions_),
      time_step_updating(time_step_updating_),
      dimension(dimension_),
      n_threads(n_threads_),
      pos(pos_) {
    assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
    assert(lower_bound < upper_bound);
    assert(beta_adjust_factor > 0.0 && beta_adjust_factor < 1.0);
    assert(num_positions > 0);

    std::vector<std::vector<double>> candidate(num_positions);
    std::vector<double> delt(num_positions);

    for(size_t i = 0; i < num_positions; i++){
        delt[i] = 0.0;
        candidate[i] = std::vector<double>(dimension);

        for(size_t j = 0; j < dimension; j++){
            candidate[i][j] = 0.0;
        }
    }
 
    candidate_positions = candidate;
    delta = delt;
}


void MultiStochasticTunnelling::iteration(const size_t seed, const size_t k){
  double sigma = compute_sigma(k);
  

  if(k % time_step_updating == 0 && k!=0){
    double best_fit = std::numeric_limits<double>::infinity();
    std::vector<double> best_location(dimension);

    #pragma omp parallel for reduction(min:best_fit)

    for(size_t i = 0; i < num_positions; i++){
        #pragma omp critical
        {
          if(pos[i].f0 < best_fit){
              best_fit = pos[i].f0;
              best_location = pos[i].best_position;
          }
        }
    }
    #pragma omp parallel for schedule(static) num_threads(n_threads) shared(best_location, best_fit)
    for(size_t i = 0; i < num_positions; i++){
        pos[i].update_best_position(best_location, best_fit);
    }
  }

  #pragma omp parallel num_threads(n_threads) shared(seed, sigma, upper_bound, lower_bound, delta, pos, candidate_positions)
  {
    const size_t thread_id = omp_get_thread_num();
    #pragma omp for schedule(static)
    for(size_t i = 0; i < num_positions; i++){


      candidate_positions[i] = pos[i].generate_new_position(lower_bound, upper_bound, seed + thread_id, sigma);
      delta[i] = mapped_function_value(candidate_positions[i], i) - mapped_function_value(pos[i].position, i);




      if(delta_condition(delta[i]) or metropolis_condition(delta[i], seed + thread_id, pos[i].beta, func(candidate_positions[i]) - func(pos[i].position), func(pos[i].best_position) - func(pos[i].position))){
          pos[i].update_window_tunnelling(1);
          pos[i].update_position(candidate_positions[i], func);
      }else{
          pos[i].update_window_tunnelling(0);
      }
      pos[i].update_betan(beta_adjust_factor, beta_thresholding);

    }    
  }
}



double MultiStochasticTunnelling::mapped_function_value(const std::vector<double>& posi, size_t i){  
  return (1.0 - std::exp(-gamma * (func(posi) - pos[i].f0)));
}

bool MultiStochasticTunnelling::delta_condition(double delt){
  if(delt <= 0){
    return true;
  }else{
    return false;
  }
}

bool MultiStochasticTunnelling::metropolis_condition(const double delta_f_stun, const size_t seed, const double beta, const double delta_f, const double old_delta) {
        static std::mt19937 gen(seed); 
        
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        double random_value = dist(gen);

        if(gamma*delta_f >= 1.e-6){

          double exp_value = std::exp(-beta * delta_f_stun);


          bool ris = random_value < exp_value;



          return ris;
        }else{
          
          double exp_value = std::exp(-beta * gamma * std::exp(gamma * old_delta));

        
          bool ris = random_value < exp_value;

          return ris;
        }
}

double MultiStochasticTunnelling::compute_sigma(size_t i){
  return sigma_max - (((sigma_max - sigma_min)/ max_iter))*i;
}


