#include <iomanip>
#include <iostream>
#include "benchmark_low_dim.hpp"

using namespace fdapde;

int main() {
	std::cout << std::setprecision(2) << std::scientific;

	LBFGS<Dynamic, WolfeLineSearch> lbfgs{500, 1e-5, 1e-2, 30};
	print_optim_benchmark(lbfgs, "LBFGS-30, Wolfe");
	
	return 0;
}