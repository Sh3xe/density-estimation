#include "benchmark_low_dim.hpp"
#include <iostream>
#include <fstream>

using namespace fdapde;

static constexpr int MAX_ITER = 500;
static constexpr double TOL = 1e-5;
static constexpr double STEP = 1e-2;

void lowdim::print_performances( const lowdim::BenchmarkResult &result, std::ostream &os ) {
	// print a markdown table, to be exported to pdf for better visibility
	os <<
		"|" << result.title << 
		"|" << result.nit_mean << 
		"|" << result.nit_std << 
		"|" << result.x_diff_mean <<
		"|" << result.x_diff_std <<
		"|" << result.f_diff_mean <<
		"|" << result.f_diff_std <<
		"|" << utils::duration_to_str(result.duration) << 
		"|" << utils::duration_to_str(result.duration_std) << 
		"|\n";

		// print(f"{title}|{nit:.4}|{nit_std:.4}|{x_diff:.2e}|{x_diff_std:.2e}|{f_diff:.2e}|{f_diff_std:.2e}|{duration_to_str(duration)}|{duration_to_str(duration_std)}")
}

void lowdim::print_optim_benchmark(
	const std::vector<lowdim::BenchmarkResult> &res,
	const std::string &title,
	std::ostream &os 
) {
	// Markdown header
	os << "### " << title << "\n";
	os << "| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std\n";
	os << "|-|-|-|-|-|-|-|-|-|\n";

	for(const auto &result: res) {
		print_performances(result, os);
	}
}

void lowdim::full_benchmark() {

	std::ofstream file{ path(ROOT_DIR) / path("outputs/cpp_bmrk_table.md"), std::ios::out | std::ios::trunc};

	if(!file) {
		std::cerr << "Cannot write to \"outputs/cpp_bmrk_table.md\"";
		return;
	}
	file << std::scientific << std::setprecision(2);
	std::cout << "Optimizing with lbfgs30" << std::endl;
	LBFGS<Dynamic> lbfgs30 {MAX_ITER, TOL, STEP, 30};
	auto lbfgs30_res = benchmark_optimizer(lbfgs30, "lbfgs30");
	print_optim_benchmark(lbfgs30_res, "lbfgs30", file);
	std::cout << "Done optimizing with lbfgs30" << std::endl;

	std::cout << "Optimizing with cg_fr" << std::endl;
	FletcherReevesCG<Dynamic> cg_fr {MAX_ITER, TOL, STEP};
	auto cg_fr_res = benchmark_optimizer(cg_fr, "cg_fr");
	print_optim_benchmark(cg_fr_res, "cg_fr", file);
	std::cout << "Done ptimizing with cg_fr" << std::endl;
	
	std::cout << "Optimizing with cg_pr" << std::endl;
	PolakRibiereCG<Dynamic> cg_pr {MAX_ITER, TOL, STEP};
	auto cg_pr_res = benchmark_optimizer(cg_pr, "cg_pr");
	print_optim_benchmark(cg_pr_res, "cg_pr", file);
	std::cout << "Done ptimizing with cg_pr" << std::endl;
	
	std::cout << "Optimizing with nelder_mead" << std::endl;
	NelderMead<Dynamic> nelder_mead {MAX_ITER, TOL};
	auto nelder_mead_res = benchmark_optimizer(nelder_mead, "nelder_mead");
	print_optim_benchmark(nelder_mead_res, "nelder_mead", file);
	std::cout << "Done ptimizing with nelder_mead" << std::endl;
}