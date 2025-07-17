#pragma once

#include <fdaPDE/optimization.h>
#include <chrono>
#include <Eigen/Dense>

#include "utilities.hpp"

class PrintValues {
public:
	PrintValues() = default;
	
	template <typename Opt, typename Obj> bool pre_update_step(Opt& opt, Obj& obj) {
		std::cout << "----- NEW ITERATION -----" << std::endl;
		std::cout << opt.x_old << std::endl;
		return false;
	}
};

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

struct LowDimBenchmarkResult {
	std::string title;
	microsec duration;
	double real_min;
	Eigen::VectorXd real_argmin;
	double opt_min;
	Eigen::VectorXd opt_argmin;
	size_t iterations;
};

void print_performances( const LowDimBenchmarkResult &result );

template<
	typename Optimizer,
	typename FieldFunc
>
LowDimBenchmarkResult benchmark_function(
	Optimizer &opt,
	FieldFunc &func,
	const Eigen::VectorXd &init_pt,
	Eigen::VectorXd real_argmin,
	double real_min,
	const std::string &name
) {
	auto begin_time = std::chrono::high_resolution_clock::now();
	
	for(int i = 0; i < 100; ++i) {
		opt.optimize(func, init_pt);
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	double duration = microsec(end_time - begin_time).count();
	
	LowDimBenchmarkResult res;
	res.title = name;
	res.duration = microsec(duration / 100.0);
	res.real_argmin = real_argmin;
	res.real_min = real_min;
	res.opt_min = opt.value();
	res.iterations = opt.n_iter();
	res.opt_argmin = opt.optimum();
	return res;
}

template <typename OptimizerType>
std::vector<LowDimBenchmarkResult> benchmark_optimizer(OptimizerType &opt, const std::string &title) {

	// Optimization functions and their init point
	fdapde::ScalarField<fdapde::Dynamic, SchafferF6>     schaffer_f6(2);
	fdapde::ScalarField<fdapde::Dynamic, Schwefel>       schwefel_2d(2);
	fdapde::ScalarField<fdapde::Dynamic, Schwefel>       schwefel_10d(10);
	fdapde::ScalarField<fdapde::Dynamic, Sphere>         sphere_2d(2);
	fdapde::ScalarField<fdapde::Dynamic, Sphere>         sphere_30d(30);
	fdapde::ScalarField<fdapde::Dynamic, Rastrigin>      rastrigin_30d(30);
	fdapde::ScalarField<fdapde::Dynamic, Rastrigin>      rastrigin_2d(2);
	fdapde::ScalarField<fdapde::Dynamic, Rosenbrock>     rosenbrock(2);
	
	// Initial points
	Eigen::VectorXd init_1 = Eigen::VectorXd::Constant(2, 1, 1.0);
	Eigen::VectorXd init_2 = Eigen::VectorXd::Constant(10, 1, 1.0);
	Eigen::VectorXd init_3 = Eigen::VectorXd::Constant(30, 1, 1.0);
	Eigen::VectorXd init_4 (2);
	init_4 << -1.2, 1;

	// Benchmark the method on different functions and print the result to markdown in the standard output
	std::vector<LowDimBenchmarkResult> res;
	res.reserve(8);

	res.push_back(benchmark_function(
		opt,
		sphere_2d,
		init_1, 
		Sphere::argmin(2),
		Sphere::min(),
		"Sphere 2D"
	));

	res.push_back(benchmark_function(
		opt,
		sphere_30d,
		init_3, 
		Sphere::argmin(30),
		Sphere::min(),
		"Sphere 30D"
	));

	res.push_back(benchmark_function(
		opt,
		schwefel_2d,
		init_1 * 20.0,
		Schwefel::argmin(2),
		Schwefel::min(),
		"Schwefel 2D" 
	));
	
	res.push_back(benchmark_function(
		opt,
		schwefel_10d,
		init_2 * 20.0,
		Schwefel::argmin(10),
		Schwefel::min(),
		"Schwefel 10D"
	));

	res.push_back(benchmark_function(
		opt,
		rastrigin_2d,
		init_1 * 3.0,
		Rastrigin::argmin(2),
		Rastrigin::min(),
		"Rastrigin 2D"
	));

	res.push_back(benchmark_function(
		opt,
		rastrigin_30d,
		init_3 * 3.0,
		Rastrigin::argmin(30),
		Rastrigin::min(),
		"Rastrigin 30D"
	));

	res.push_back(benchmark_function(
		opt,
		schaffer_f6,
		init_1 * 5.0,
		SchafferF6::argmin(2),
		SchafferF6::min(),
		"Schaffer F6"
	));

	res.push_back(benchmark_function(
		opt,
		rosenbrock,
		init_4, 
		Rosenbrock::argmin(2),
		Rosenbrock::min(),
		"Rosenbrock 2D"
	));


	return res;
}

template <typename OptimizerType>
void print_optim_benchmark(OptimizerType &opt, const std::string &title) {
	// Markdown header
	std::cout << "### `" << title << "`\n";
	std::cout << "|Method|n_iter| x-x* l2 err | f-f* l2 err | time |" << std::endl;
	std::cout << "|-|-|-|-|-|" << std::endl;

	auto benchmark = benchmark_optimizer(opt, title);

	for(const auto &result: benchmark) {
		print_performances(result);
	}
}