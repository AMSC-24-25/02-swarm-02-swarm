#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <mpi.h>

#include "DistributedMultiStochasticTunnelling.hpp"


DistributedMultiStochasticTunnelling::DistributedMultiStochasticTunnelling(const int world_rank_, const int world_size_,std::vector<Position>& pos_, const double lower_bound_, const double upper_bound_,const double sigma_max_, const double sigma_min_,
                        const double gamma_, const double beta_adjust_factor_,
                        const size_t max_iter_, const ObjectiveFunction& func_, double beta_thresholding_, const size_t num_positions_, const size_t time_step_updating_, const size_t dimension_)
    : 
      world_rank(world_rank_),
      world_size(world_size_),
      lower_bound(lower_bound_),
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
      pos(pos_) {
	assert(world_size > 0);
	assert(world_rank >= 0 && world_rank <= world_size - 1);
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


void DistributedMultiStochasticTunnelling::iteration(const size_t seed, const size_t k){
  double sigma = compute_sigma(k);
    const size_t local_size = num_positions / world_size;
    const size_t remainder = num_positions % world_size;
    const size_t start = (static_cast<size_t>(world_rank) < remainder) ? 
                        world_rank * (local_size + 1) : 
                        remainder * (local_size + 1) + (world_rank - remainder) * local_size;
    const size_t end = start + local_size + (static_cast<size_t>(world_rank) < remainder ? 1 : 0);

    MPI_Barrier(MPI_COMM_WORLD);

    if(k % time_step_updating == 0 && k != 0) {
        double local_best_fit = std::numeric_limits<double>::infinity();
        std::vector<double> local_best_location(dimension);
        
        for(size_t i = start; i < end; i++) {
            if(pos[i].f0 < local_best_fit) {
                local_best_fit = pos[i].f0;
                local_best_location = pos[i].best_position;
            }
        }

        if(world_rank == 0) {
            std::vector<double> all_fits(world_size);
            std::vector<std::vector<double>> all_locations(world_size, std::vector<double>(dimension));
            
            for(int i = 1; i < world_size; i++) {
                MPI_Recv(&all_fits[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(all_locations[i].data(), dimension, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            all_fits[0] = local_best_fit;
            all_locations[0] = local_best_location;

            auto min_it = std::min_element(all_fits.begin(), all_fits.end());
            int best_rank = std::distance(all_fits.begin(), min_it);
            std::vector<double> global_best = all_locations[best_rank];
            double global_best_fit = *min_it;

            for(int i = 1; i < world_size; i++) {
                MPI_Send(&global_best_fit, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(global_best.data(), dimension, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }

            for(size_t i = start; i < end; i++) {
                pos[i].update_best_position(global_best, global_best_fit);
            }

        } else {
            MPI_Send(&local_best_fit, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(local_best_location.data(), dimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            double global_best_fit;
            std::vector<double> global_best(dimension);
            MPI_Recv(&global_best_fit, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(global_best.data(), dimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(size_t i = start; i < end; i++) {
                pos[i].update_best_position(global_best, global_best_fit);
            }
        }
    }



    std::mt19937 gen(seed + world_rank);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for(size_t i = start; i < end; i++) {
        candidate_positions[i] = pos[i].generate_new_position(
            lower_bound, upper_bound, 
            gen(), 
            sigma
        );

        delta[i] = mapped_function_value(candidate_positions[i], i) - 
                  mapped_function_value(pos[i].position, i);

        if(delta_condition(delta[i], i) || 
           metropolis_condition(delta[i], gen(), pos[i].beta,
                              func(candidate_positions[i]) - func(pos[i].position),
                              func(pos[i].best_position) - func(pos[i].position), i)) 
        {
            pos[i].update_position(candidate_positions[i], func);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}



double DistributedMultiStochasticTunnelling::mapped_function_value(const std::vector<double>& posi, size_t i){  
  return (1.0 - std::exp(-gamma * (func(posi) - pos[i].f0)));
}

bool DistributedMultiStochasticTunnelling::delta_condition(double delt, size_t i){
  if(delt <= 0){
    pos[i].update_window_tunnelling(1);
    pos[i].update_betan(beta_adjust_factor, beta_thresholding);
    return true;
  }else{
    pos[i].update_window_tunnelling(0);
    pos[i].update_betan(beta_adjust_factor, beta_thresholding);
    return false;
  }
}

bool DistributedMultiStochasticTunnelling::metropolis_condition(const double delta_f_stun, const size_t seed, const double beta, const double delta_f, const double old_delta, size_t i) {
        static std::mt19937 gen(seed); 
        
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        double random_value = dist(gen);

        //std::cout<<"gamma * delta_f: "<<gamma*delta_f<<std::endl;

        if(gamma*delta_f >= 1.e-6){

          double exp_value = std::exp(-beta * delta_f_stun);

          //std::cout<<"beta: "<<beta<<std::endl;

          //std::cout<<"percentage: "<<exp_value<<std::endl;

          bool ris = random_value < exp_value;

          if(ris == true){
            pos[i].update_window_tunnelling(1);
          }else{
            pos[i].update_window_tunnelling(0);
          }

          pos[i].update_betan(beta_adjust_factor, beta_thresholding);

          return ris;
        }else{
          
          double exp_value = std::exp(-beta * gamma * std::exp(gamma * old_delta));

          //std::cout<<"percentage approximated: "<<exp_value<<std::endl;
        
          bool ris = random_value < exp_value;

          if(ris == true){
            pos[i].update_window_tunnelling(1);
          }else{
            pos[i].update_window_tunnelling(0);
          }

          pos[i].update_betan(beta_adjust_factor, beta_thresholding);

          return ris;
        }
}

double DistributedMultiStochasticTunnelling::compute_sigma(size_t i){
  return sigma_max - (((sigma_max - sigma_min)/ max_iter))*i;
}

void DistributedMultiStochasticTunnelling::update_beta_thresholding(size_t k){
  if(k == 1){
    beta_thresholding = beta_thresholding * 2/3;
  }

  if(k==2){
    beta_thresholding = beta_thresholding / 2;
  }
}