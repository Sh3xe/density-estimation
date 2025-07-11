#pragma once

#include "utilities.hpp"

struct DEBenchmarkResult {
	size_t iterations;
	vector_t log_density;
	double loss;
	microsec duration;
	std::string title;
};

template< typename Optimizer >
DEBenchmarkResult benchmark_solver(
	const fdapde::DEPDE &solver,
	const Optimizer &optimizer,
	const vector_t &g_init,
	const std::string &title )
{
	auto begin_time = std::chrono::high_resolution_clock::now();
	
	auto end_time = std::chrono::high_resolution_clock::now();
	
	DEBenchmarkResult res;
	res.title = title;
	res.duration = microsec(end_time - begin_time);
	res.iterations = opt.n_iter();
	return res;
}

std::vector<DEBenchmarkResult> benchmark_de_problem(
	const std::string &mesh_dir_str,
	const std::string &title
);