#include "benchmark_low_dim.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ostream>
#include <fdaPDE/optimization.h>
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include "utilities.hpp"

using namespace fdapde;
using namespace Eigen;

static constexpr int MAX_ITER = 500;
static constexpr double TOL = 1e-5;
static constexpr double STEP = 1e-2;

struct SchafferF6 {
	SchafferF6() = default;

	static VectorXd argmin(int dim) {
		return VectorXd::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const VectorXd& x) const {
		double r = x.dot(x);
		double a = std::sin(std::sqrt(r));
		double b = 1.0 + 0.001 * r;
		return 0.5 + ( a*a - 0.5 ) / ( b*b );
	}
};

struct Schwefel {
	Schwefel() = default;

	static VectorXd argmin(int dim) {
		return VectorXd::Constant(dim, 1, 420.9687);
	}

	static double min() { return 0.0; }

	double operator()(const VectorXd& x) const {
		double s = 0.0;
		for(int i = 0; i < x.rows(); ++i) {
			if(x[i] > 500 || x[i] < -500)
				s += 0.02 * x[i] * x[i];
			else
				s += -x[i] * std::sin(std::sqrt(std::abs(x[i])));
		}
		double dimension = (double)x.rows();
		return 418.9829 * dimension + s;
	}
};

struct Rosenbrock {
	Rosenbrock() = default;

	static VectorXd argmin(int dim) {
		return VectorXd::Constant(dim, 1, 1.0);
	}

	static double min() { return 0.0; }

	double operator()(const VectorXd& x) const {
		constexpr double a = 1.0;
		constexpr double b = 100.0;
		return (a-x[0])*(a-x[0]) + b*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
	}
};

struct Rastrigin {
	Rastrigin() = default;

	static VectorXd argmin(int dim) {
		return VectorXd::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const VectorXd& x) const {
		double s = 0.0;
		for(int i = 0; i < x.rows(); ++i) {
			if(x[i] > 5.12 || x[i] < -5.12)
				s += 10.0 * x[i] * x[i];
			else
				s += x[i] * x[i] - 10.0*std::cos(2*3.14159*x[i]);
		}
		double dimension = (double)x.rows();
		return dimension * 10 + s;
	}
};

struct Sphere {
	Sphere() = default;

	static VectorXd argmin(int dim) {
		return VectorXd::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const VectorXd& x) const {
		return x.dot(x);
	}
};

struct EvalCounterHook {
	EvalCounterHook() {reset();}

	void reset() {
		obj_eval_count = 0;
		grad_eval_count = 0;
	}

	template <typename Opt, typename Obj> bool exec_grad_hooks(Opt& opt, Obj& obj) {
		++grad_eval_count;
		return false;
	}

	template <typename Opt, typename Obj> bool exec_eval_hooks(Opt& opt, Obj& obj) {
		++obj_eval_count;
		return false;
	}

	int obj_eval_count {0};
	int grad_eval_count {0};
};

struct GeneticSave {
	GeneticSave() {}

	int current = 0;
	std::vector<std::string> columns;

	void reset(int n) {
	}

	template <typename Opt> bool select_hook(Opt& opt) {
		++current;
		if(columns.size() == 0) reset(opt.population.cols());
		utils::write_csv( "genetic_csv/" + std::to_string(current) + std::string("_pop.csv"), {"x", "y"}, opt.population.transpose() );
		return false;
	}

	template <typename Opt, typename Obj> bool eval_hook(Opt& opt, Obj& obj) {
		return false;
	}
};

struct LDBenchmarkResult {
	std::string title;
	microsec duration;
	microsec duration_std;
	double x_diff_mean;
	double x_diff_std;
	double f_diff_mean;
	double f_diff_std;
	double nit_mean;
	double nit_std;
};

template<
	typename Optimizer,
	typename FieldFunc,
	typename ...Hooks
>
LDBenchmarkResult benchmark_function(
	Optimizer &opt,
	FieldFunc &func,
	const MatrixXd &init_pts,
	VectorXd real_argmin,
	double real_min,
	const std::string &title,
	const std::string &opt_title,
	bool output_csv,
	Hooks &&...hooks
) {
	// [nit, x_diff, f_diff, p_id, duration_microsec]
	MatrixXd df = MatrixXd::Zero(init_pts.rows(), 7);
	EvalCounterHook eval_count_hook;

	for(int i = 0; i < init_pts.rows(); ++i) {
		// Perform benchmark on one initial point
		eval_count_hook.reset();
		auto begin_time = std::chrono::high_resolution_clock::now();
		opt.optimize(func, init_pts.row(i), eval_count_hook, std::forward<Hooks>(hooks)...);
		auto end_time = std::chrono::high_resolution_clock::now();
		double duration = microsec(end_time - begin_time).count();

		// Save the results to csv
		df(i, 0) = (double)opt.n_iter();
		df(i, 1) = (real_argmin - opt.optimum()).norm();
		df(i, 2) = std::abs(real_min - opt.value());
		df(i, 3) = (double)i;
		df(i, 4) = duration;
		df(i, 5) = (double)eval_count_hook.obj_eval_count;
		df(i, 6) = (double)eval_count_hook.grad_eval_count;
	}

	if(output_csv) {
		std::string file_title = "outputs/cpp_bmrk_" + opt_title + "_" + title + "_" + std::to_string(init_pts.rows()) + "d.csv";
		
		utils::write_csv(
			file_title,
			{"nit","x_diff","f_diff","p_id","duration_microsec", "obj_eval", "grad_eval"},
			df
		);
	}

	LDBenchmarkResult res;
	res.title = title;
	res.nit_mean = utils::mean(df, 0);
	res.nit_std = utils::std(df, 0);
	res.x_diff_mean = utils::mean(df, 1);
	res.x_diff_std = utils::std(df, 1);
	res.f_diff_mean = utils::mean(df, 2);
	res.f_diff_std = utils::std(df, 2);
	res.duration = microsec(utils::mean(df, 4));
	res.duration_std = microsec(utils::std(df, 4));

	return res;
}

template <typename OptimizerType, typename ...Hooks>
std::vector<LDBenchmarkResult> benchmark_optimizer(
	OptimizerType &opt,
	const std::string &title,
	bool output_csv,
	Hooks &&...hooks
) {
	std::cout << "Benchmarking " << title << std::endl;
	
	// Benchmark the method on different functions and print the result to markdown in the standard output
	std::vector<int> dims {2,10,30};
	std::vector<LDBenchmarkResult> res; res.reserve(11);
	 
	for(auto &d: dims) {
		fdapde::ScalarField<fdapde::Dynamic, Sphere> shpere_func(d);
		res.push_back(benchmark_function(
			opt, shpere_func,
			utils::load_csv("data/lowdim_inits/" + std::to_string(d) + "_sphere.csv"), 
			Sphere::argmin(d), Sphere::min(), "sphere_" + std::to_string(d) + "d", title,
			output_csv, std::forward<Hooks>(hooks)...
		));
	}

	for(auto &d: dims) {
		fdapde::ScalarField<fdapde::Dynamic, Schwefel> schwefel_func(d);
		res.push_back(benchmark_function(
			opt, schwefel_func,
			utils::load_csv("data/lowdim_inits/" + std::to_string(d) + "_schwefel.csv"),
			Schwefel::argmin(d), Schwefel::min(), "schwefel_" + std::to_string(d) + "d", title,
			output_csv, std::forward<Hooks>(hooks)...
		));
	}

	for(auto &d: dims) {
		fdapde::ScalarField<fdapde::Dynamic, Rastrigin> rastrigin_func(d);
		res.push_back(benchmark_function(
			opt, rastrigin_func,
			utils::load_csv("data/lowdim_inits/" + std::to_string(d) + "_rastrigin.csv"),
			Rastrigin::argmin(d), Rastrigin::min(), "rastrigin_" + std::to_string(d) + "d", title,
			output_csv, std::forward<Hooks>(hooks)...
		));
	}

	fdapde::ScalarField<fdapde::Dynamic, SchafferF6> schaffer_f6(2);
	res.push_back(benchmark_function(
		opt, schaffer_f6,
		utils::load_csv("data/lowdim_inits/schaffer_f6.csv"),
		SchafferF6::argmin(2), SchafferF6::min(), "schaffer_f6", title,
		output_csv, std::forward<Hooks>(hooks)...
	));
	
	fdapde::ScalarField<fdapde::Dynamic, Rosenbrock> rosenbrock(2);
	res.push_back(benchmark_function(
		opt, rosenbrock,
		utils::load_csv("data/lowdim_inits/rosenbrock.csv"), 
		Rosenbrock::argmin(2),Rosenbrock::min(), "rosenbrock", title,
		output_csv, std::forward<Hooks>(hooks)...
	));

	return res;
}

void print_performances( const LDBenchmarkResult &result, std::ostream &os ) {
	// print a markdown table, to be exported to pdf for better visibility
	os <<
		"|" << result.title << std::fixed <<
		"|" << result.nit_mean << 
		"|" << result.nit_std << std::scientific <<
		"|" << result.x_diff_mean << 
		"|" << result.x_diff_std <<
		"|" << result.f_diff_mean <<
		"|" << result.f_diff_std << std::fixed <<
		"|" << utils::duration_to_str(result.duration) << 
		"|" << utils::duration_to_str(result.duration_std) << 
		"|\n";
}

void print_optim_benchmark(
	const std::vector<LDBenchmarkResult> &res,
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

void lowdim_full_benchmark(bool output_csv) {

	std::ofstream file{ path(ROOT_DIR) / path("outputs/cpp_bmrk_table.md"), std::ios::out | std::ios::trunc};

	if(!file) {
		std::cerr << "Cannot write to \"outputs/cpp_bmrk_table.md\"";
		return;
	}

	GeneticOptim<fdapde::Dynamic> go {100, 1e-3, 5, 30};
	fdapde::ScalarField<fdapde::Dynamic, Schwefel> schwefel_func(2);
	MatrixXd init_points = utils::load_csv("data/lowdim_inits/2_rastrigin.csv");
	go.optimize(schwefel_func, init_points.row(1), BinaryTournamentSelection(), GaussianMutation(4.0, 0.95), GeneticSave());

	// {
	// 	BFGS<fdapde::Dynamic> bfgs {MAX_ITER, TOL, STEP};
	// 	auto bfgs_res = benchmark_optimizer(bfgs, "bfgs", output_csv, WolfeLineSearch());
	// 	print_optim_benchmark(bfgs_res, "bfgs", file);
	// }

	// {
	// 	LBFGS<fdapde::Dynamic> lbfgs30 {MAX_ITER, TOL, STEP, 30};
	// 	auto lbfgs30_res = benchmark_optimizer(lbfgs30, "lbfgs30", output_csv, WolfeLineSearch());
	// 	print_optim_benchmark(lbfgs30_res, "lbfgs30", file);
	// }

	// {
	// 	FletcherReevesCG<fdapde::Dynamic> cg_fr {MAX_ITER, TOL, STEP};
	// 	auto cg_fr_res = benchmark_optimizer(cg_fr, "cg_fr", output_csv, WolfeLineSearch());
	// 	print_optim_benchmark(cg_fr_res, "cg_fr", file);
	// }
		
	// {
	// 	PolakRibiereCG<fdapde::Dynamic> cg_pr {MAX_ITER, TOL, STEP, true};
	// 	auto cg_pr_res = benchmark_optimizer(cg_pr, "cg_pr_restart", output_csv, WolfeLineSearch());
	// 	print_optim_benchmark(cg_pr_res, "cg_pr_restart", file);
	// }

	// {
	// 	PolakRibierePlsCG<fdapde::Dynamic> cg_prp {MAX_ITER, TOL, STEP, true};
	// 	auto cg_prp_res = benchmark_optimizer(cg_prp, "cg_prp_restart", output_csv, WolfeLineSearch());
	// 	print_optim_benchmark(cg_prp_res, "cg_prp_restart", file);
	// }

	// {
	// 	NelderMead<fdapde::Dynamic> nelder_mead {MAX_ITER, TOL};
	// 	auto nelder_mead_res = benchmark_optimizer(nelder_mead, "nelder_mead", output_csv);
	// 	print_optim_benchmark(nelder_mead_res, "nelder_mead", file);
	// }

	// {
	// 	GeneticOptim<fdapde::Dynamic> genetic_bin_co {MAX_ITER, TOL, 5, 30};
	// 	auto genetic_bin_co_res = benchmark_optimizer(genetic_bin_co, "genetic_bin_co", output_csv, BinaryTournamentSelection(),GaussianMutation());
	// 	print_optim_benchmark(genetic_bin_co_res, "genetic_bin_co", file);
	// }

	// {
	// 	GeneticOptim<fdapde::Dynamic> genetic_rk_gaus {MAX_ITER, TOL, 5, 30};
	// 	auto genetic_rk_gaus_res = benchmark_optimizer(genetic_rk_gaus, "genetic_rk_gaus", output_csv, RankSelection(), GaussianMutation());
	// 	print_optim_benchmark(genetic_rk_gaus_res, "genetic_rk_gaus", file);
	// }
}