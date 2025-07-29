#pragma once

#include <fdaPDE/optimization.h>
#include <Eigen/Dense>
#include "utilities.hpp"

namespace lowdim {

// Test functions

struct SchafferF6 {
	SchafferF6() = default;

	static Eigen::VectorXd argmin(int dim) {
		return Eigen::VectorXd::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const Eigen::VectorXd& x) const {
		double r = x.dot(x);
		double a = std::sin(std::sqrt(r));
		double b = 1.0 + 0.001 * r;
		return 0.5 + ( a*a - 0.5 ) / ( b*b );
	}
};

struct Schwefel {
	Schwefel() = default;

	static Eigen::VectorXd argmin(int dim) {
		return Eigen::VectorXd::Constant(dim, 1, 420.9687);
	}

	static double min() { return 0.0; }

	double operator()(const Eigen::VectorXd& x) const {
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

	static Eigen::VectorXd argmin(int dim) {
		return Eigen::VectorXd::Constant(dim, 1, 1.0);
	}

	static double min() { return 0.0; }

	double operator()(const Eigen::VectorXd& x) const {
		constexpr double a = 1.0;
		constexpr double b = 100.0;
		return (a-x[0])*(a-x[0]) + b*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
	}
};

struct Rastrigin {
	Rastrigin() = default;

	static Eigen::VectorXd argmin(int dim) {
		return Eigen::VectorXd::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const Eigen::VectorXd& x) const {
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

	static Eigen::VectorXd argmin(int dim) {
		return Eigen::VectorXd::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const Eigen::VectorXd& x) const {
		return x.dot(x);
	}
};

// Benchmark result

struct BenchmarkResult {
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

// Utils (printing &optimizing)

void print_optim_benchmark(const std::vector<BenchmarkResult> &res, const std::string &title, std::ostream &os);

void print_performances( const BenchmarkResult &result, std::ostream &os );

void full_benchmark(bool output_csv = false);

// Optimization benchmark

template<
	typename Optimizer,
	typename FieldFunc,
	typename ...Hooks
>
BenchmarkResult benchmark_function(
	Optimizer &opt,
	FieldFunc &func,
	const Eigen::MatrixXd &init_pts,
	Eigen::VectorXd real_argmin,
	double real_min,
	const std::string &title,
	const std::string &opt_title,
	bool output_csv,
	Hooks &&...hooks
) {
	// [nit, x_diff, f_diff, p_id, duration_microsec]
	Eigen::MatrixXd df = Eigen::MatrixXd::Zero(init_pts.rows(), 5);

	for(int i = 0; i < init_pts.rows(); ++i) {
		// Perform benchmark on one initial point
		auto begin_time = std::chrono::high_resolution_clock::now();
		opt.optimize(func, init_pts.row(i), std::forward<Hooks>(hooks)...);
		auto end_time = std::chrono::high_resolution_clock::now();
		double duration = microsec(end_time - begin_time).count();

		// Save the results to csv
		df(i, 0) = (double)opt.n_iter();
		df(i, 1) = (real_argmin - opt.optimum()).norm();
		df(i, 2) = std::abs(real_min - opt.value());
		df(i, 3) = (double)i;
		df(i, 4) = duration;
	}

	if(output_csv) {
		std::string file_title = "outputs/cpp_bmrk_" + opt_title + "_" + title + "_" + std::to_string(init_pts.rows()) + "d.csv";
		
		utils::write_csv(
			file_title,
			{"nit","x_diff","f_diff","p_id","duration_microsec"},
			df
		);
	}

	BenchmarkResult res;
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
std::vector<BenchmarkResult> benchmark_optimizer(
	OptimizerType &opt,
	const std::string &title,
	bool output_csv,
	Hooks &&...hooks
) {
	std::cout << "Benchmarking " << title << std::endl;
	
	// Benchmark the method on different functions and print the result to markdown in the standard output
	std::vector<int> dims {2,10,30};
	std::vector<BenchmarkResult> res; res.reserve(11);

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

} // namespace lowdim