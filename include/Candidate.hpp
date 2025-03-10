#ifndef CANDIDATE_HPP
#define CANDIDATE_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <vector>

#include "ObjectiveFunction.hpp"

class Candidate {
	public:
	double f0;
	std::vector<double> candidate;

	private:
	const ObjectiveFunction& func;

	public:
	Candidate(const size_t dimension, const double lower_bound, const double upper_bound, const size_t seed, const ObjectiveFunction& func);
	Candidate(const std::vector<double> candidate_, const ObjectiveFunction& func);
	void updatePosition(const std::vector<double>& position);

};

#endif