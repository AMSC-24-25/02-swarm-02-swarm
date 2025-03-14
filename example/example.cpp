#include <iostream>
#include <vector>
#include <cassert>
#include <numeric>

#include <ObjectiveFunction.hpp>
#include <Algorithm.hpp>

class MyCustomFunction : public ObjectiveFunction {
   public:
	double operator()(const std::vector<double>& position) const {
		assert(position.size() > 0);

		return std::reduce(position.begin(), position.end());
	}
};

int main(void) {
	std::cout << "Minimizing a custom function" << std::endl;
	std::cout << std::endl;

	const size_t dimensions = 2;
	const size_t num_creatures = 1'000;
	const size_t max_iterations = 100;
	const size_t seed = 42;
	const double lower_bound = -10.0;
	const double upper_bound = 10.0;
	const size_t n_threads = 1;
	const double mutation_rate = 0.2;
	const double survival_rate = 0.5;

	std::unique_ptr<ObjectiveFunction> mcf = std::make_unique<MyCustomFunction>();

	const std::pair<std::vector<double>, double> result =
		algorithm::run_genetic_openmp(dimensions, num_creatures, max_iterations, seed, lower_bound, upper_bound,
									  mutation_rate, survival_rate, mcf, n_threads, true);

	std::cout << "Minimum found:" << std::endl;
	std::cout << "  f(";
	for (size_t i{0}; i < result.first.size(); i++) {
		std::cout << result.first.at(i);
		if (i < result.first.size() - 1) {
			std::cout << ", ";
		}
	}
	std::cout << ") = " << result.second << std::endl;
	std::cout << std::endl;

	return 0;
}