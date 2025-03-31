#include <vector>
#include <limits>
#include <random>
#include <algorithm>
#include <cassert>

#include "Position.hpp"

Position::Position(const size_t dimensions_, const double lower_bound, const double upper_bound, const size_t seed, const double beta_, const ObjectiveFunction& func, const size_t window_tunnelling_)
    : f0(std::numeric_limits<double>::infinity()), dimensions(dimensions_), window_tunnelling(window_tunnelling_) {
    assert(dimensions_ > 0);
	assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
	assert(lower_bound < upper_bound);

    std::vector<double> position_(dimensions);

    std::mt19937 rnd{seed};
	std::uniform_real_distribution<double> position_dist{lower_bound, upper_bound};


    for (size_t i = 0; i < dimensions; i++) {
		position_[i] = position_dist(rnd);
	}

    position = position_;
    best_position = position_;
    beta = beta_;
    f0 = func(position);
    

    std::vector<int> tunnelling(window_tunnelling);
    for(size_t i = 0; i< window_tunnelling; i++){
        tunnelling[i] = 0;   
    }

    avg_tunnelling = tunnelling;
}


const std::vector<double> Position::generate_new_position(const double lower_bound,
				const double upper_bound, const size_t seed, const double sigma){

    static std::mt19937 gen(seed); 
    std::normal_distribution<double> dist(0.0, sigma);
    std::vector<double> posit(dimensions);
    
    for (size_t i = 0; i< dimensions; i++){
        posit[i] = position[i];
        double delta = dist(gen);
        posit[i] += delta;
        posit[i] = std::clamp(posit[i], lower_bound, upper_bound);
    }
    
    return posit;
}

void Position::update_position(std::vector<double> position_, const ObjectiveFunction& func){
    if(func(position_) < f0){
        f0 = func(position_);
        best_position = position_;
    }
    
    position = position_;
}


void Position::update_betan(double beta_adjust_factor, double beta_thresholding){
    std::cout<<"average tunnelling: "<<compute_avg_tunnelling() <<std::endl;
    if(compute_avg_tunnelling() > beta_thresholding) {
        beta /= beta_adjust_factor;
    }else{
        beta *= beta_adjust_factor;
    }

}

double Position::compute_avg_tunnelling(){
    double avg = 0;
    for(size_t i = 0; i < window_tunnelling; i++){
        avg += avg_tunnelling[i];
    }

    avg /= window_tunnelling;
    return avg;
}
 
void Position::update_window_tunnelling(int tunnelled){

    std::cout<<"tunneled: "<<tunnelled<<std::endl;

    for(size_t i = 0; i < window_tunnelling - 1; i++){
        avg_tunnelling[i] = avg_tunnelling[i + 1];
    }

    avg_tunnelling[window_tunnelling - 1] = tunnelled;
}

void Position::update_best_position(std::vector<double> forced_best_position, double best_fit){
    best_position = forced_best_position;
    f0 = best_fit;
}