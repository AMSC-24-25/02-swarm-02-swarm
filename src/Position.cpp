#include <vector>
#include <limits>
#include <random>
#include <algorithm>
#include <cassert>

#include "Position.hpp"

Position::Position(const size_t dimensions_, const double lower_bound, const double upper_bound, const size_t seed, const double beta_, const ObjectiveFunction& func, const int moving_avg_window_)
    : f0(std::numeric_limits<double>::infinity()), dimensions(dimensions_), moving_avg_window(moving_avg_window_) {
    assert(dimensions_ > 0);
	assert(std::isfinite(lower_bound));
	assert(std::isfinite(upper_bound));
	assert(lower_bound < upper_bound);
    assert(moving_avg_window > 0);

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
    
    std::vector<double> mean(moving_avg_window);
    for(size_t i = 0; i< moving_avg_window; i++){
        mean[i] = 0;
    }

    avg_function = mean;
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

void Position::increase_avg_window_at_position(double stun_func_value, const size_t index){
    avg_function[index] = stun_func_value;
}

void Position::update_avg_window(double stun_func_value){
    for(size_t i = 0; i < moving_avg_window - 1; i++){
        avg_function[i] = avg_function[i + 1];
    }

    avg_function[moving_avg_window - 1] = stun_func_value;
}

double Position::compute_avg_window_value(){
    double avg = 0.0;
    
    for(size_t i = 0; i < moving_avg_window; i++){
        avg += avg_function[i];
    }

    return avg /= moving_avg_window;
}

void Position::update_beta(double thresholing, double beta_adjust_factor){
    std::cout<<"the func median values is: "<<compute_avg_window_value()<<std::endl;
    if(compute_avg_window_value() > thresholing) {
        beta *= beta_adjust_factor;
        std::cout<<"beta: "<<beta<<std::endl;
    }else{
        beta /= beta_adjust_factor;
        std::cout<<"beta: "<<beta<<std::endl;
    }
}