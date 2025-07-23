#include "benchmark_low_dim.hpp"
using namespace lowdim;

void print_performances( const LowDimBenchmarkResult &result ) {
	// print a markdown table, to be exported to pdf for better visibility
	double x_err = (result.real_argmin - result.opt_argmin).norm();
	double f_err = std::abs(result.real_min - result.opt_min);

	std::cout << 
		"|" << result.title << 
		"|" << result.iterations << 
		"|" << x_err <<
		"|" << f_err <<
		"|" << duration_to_str(result.duration) << 
		"|\n";
}