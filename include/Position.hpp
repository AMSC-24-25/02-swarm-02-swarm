#ifndef POSITION_HPP
#define POSITION_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>

#include "ObjectiveFunction.hpp"

class Position {
    public:
        std::vector<double> position;
        double f0;
        const size_t dimensions;
        double beta;
        const size_t moving_avg_window;
        std::vector<double> avg_function;
        std::vector<double> best_position;
        std::vector<int> avg_tunnelling;

        const size_t window_tunnelling;

        Position(const size_t dimensions, const double lower_bound, const double upper_bound, const size_t seed, double beta, const ObjectiveFunction& func, const int moving_avg_window, const size_t window_tunnelling);

        const std::vector<double> generate_new_position(const double lower_bound,
				const double upper_bound, const size_t seed, const double sigma);

        void update_position(std::vector<double> position_, const ObjectiveFunction& func);

        void increase_avg_window_at_position(double stun_func_value, const size_t index);

        void update_avg_window(double stun_func_value);

        double compute_avg_window_value();

        void update_beta(double beta_adjust_factor);

        double compute_thresholding();

        double compute_avg_tunnelling();

        void update_betan(double beta_adjust_factor, double beta_thresholding);

        void update_window_tunnelling(int tunnelled);

    
};




#endif