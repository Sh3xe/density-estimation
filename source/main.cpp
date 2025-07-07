#include <fdaPDE/optimization.h>
#include <iomanip>
#include <chrono>
#include <Eigen/Dense>

using namespace fdapde;
using vector_t = Eigen::Matrix<double, Dynamic, 1>;
using microsec = std::chrono::duration<double, std::micro>;

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

	static vector_t argmin(int dim) {
		return vector_t::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const vector_t& x) const {
		double r = x.dot(x);
		double a = std::sin(std::sqrt(r));
		double b = 1.0 + 0.001 * r;
		return 0.5 + ( a*a - 0.5 ) / ( b*b );
	}
};

struct Schwefel {
	Schwefel() = default;

	static vector_t argmin(int dim) {
		return vector_t::Constant(dim, 1, 420.9687);
	}

	static double min() { return 0.0; }

	double operator()(const vector_t& x) const {
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

	static vector_t argmin(int dim) {
		fdapde_assert(dim == 2);
		return vector_t::Constant(dim, 1, 1.0);
	}

	static double min() { return 0.0; }

	double operator()(const vector_t& x) const {
		constexpr double a = 1.0;
		constexpr double b = 100.0;
		return (a-x[0])*(a-x[0]) + b*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
	}
};

struct Rastrigin {
	Rastrigin() = default;

	static vector_t argmin(int dim) {
		return vector_t::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const vector_t& x) const {
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

	static vector_t argmin(int dim) {
		return vector_t::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const vector_t& x) const {
		return x.dot(x);
	}
};

std::string duration_to_str(const microsec &duration) {
	if(duration.count() < 1e3) {
		return std::to_string((int)duration.count()) + std::string("μs");
	} else if(duration.count() < 1e6) {
		return std::to_string((int)(duration.count()*1e-3) ) + std::string("ms");
	} else if(duration.count() < 1e9) {
		return std::to_string((int)(duration.count()*1e-6) ) + std::string("s");
	} else {
		return ">1000s";
	}
}

template <typename OptimizerType>
void print_performances(
	OptimizerType &opt,
	const vector_t &argmin,
	double min,
	const microsec &duration,
	const std::string &test_case
) {
	// print a markdown table, to be exported to pdf for better visibility
	double x_err = (argmin - opt.optimum()).norm();
	double f_err = std::abs(min - opt.value());

	std::cout << 
		"|" << test_case << 
		"|" << opt.n_iter() << 
		"|" << x_err <<
		"|" << f_err <<
		"|" << duration_to_str(duration) << 
		"|\n";
}

template<
	typename OptimizerFunc,
	typename Optimizer,
	typename FieldFunc
>
void print_benchmark(
	Optimizer &opt,
	FieldFunc &func,
	const vector_t &init_pt,
	int dimension,
	const std::string &name
) {
	auto begin_time = std::chrono::high_resolution_clock::now();
	
	for(int i = 0; i < 100; ++i) {
		opt.optimize(func, init_pt);
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double duration = microsec(end_time - begin_time).count();
	vector_t real_argmin = OptimizerFunc::argmin(dimension);
	double real_min = OptimizerFunc::min();
	
	print_performances(
		opt,
		real_argmin,
		real_min,
		microsec(duration / 100.0),
		name
	);
}

template <typename OptimizerType>
void test_optim(OptimizerType &opt, const std::string &title) {

	// Optimization functions and their init point
	ScalarField<Dynamic, SchafferF6>     schaffer_f6(2);
	ScalarField<Dynamic, Schwefel>       schwefel_2d(2);
	ScalarField<Dynamic, Schwefel>       schwefel_10d(10);
	ScalarField<Dynamic, Sphere>         sphere_2d(2);
	ScalarField<Dynamic, Sphere>         sphere_30d(30);
	ScalarField<Dynamic, Rastrigin>      rastrigin_30d(30);
	ScalarField<Dynamic, Rastrigin>      rastrigin_2d(2);
	ScalarField<Dynamic, Rosenbrock>     rosenbrock(2);
	
	// Initial points
	vector_t init_1 = vector_t::Constant(2, 1, 1.0);
	vector_t init_2 = vector_t::Constant(10, 1, 1.0);
	vector_t init_3 = vector_t::Constant(30, 1, 1.0);
	vector_t init_4 (2);
	init_4 << -1.2, 1;

	// Markdown header
	std::cout << "### `" << title << "`\n";
	std::cout << "|Method|n_iter| x-x* l2 err | f-f* l2 err | time |" << std::endl;
	std::cout << "|-|-|-|-|-|" << std::endl;

	// Benchmark the method on different functions and print the result to markdown in the standard output
	print_benchmark<Sphere>(opt, sphere_2d, init_1, 2, "Sphere 2D" );
	print_benchmark<Sphere>(opt, sphere_30d, init_3, 30, "Sphere 30D" );
	print_benchmark<Schwefel>(opt, schwefel_2d, init_1 * 20.0, 2, "Schwefel 2D" );
	print_benchmark<Schwefel>(opt, schwefel_10d, init_2 * 20.0, 10, "Schwefel 10D" );
	print_benchmark<Rastrigin>(opt, rastrigin_2d, init_1 * 3.0, 2, "Rastrigin 2D" );
	print_benchmark<Rastrigin>(opt, rastrigin_30d, init_3 * 3.0, 30, "Rastrigin 30D" );
	print_benchmark<SchafferF6>(opt, schaffer_f6, init_1 * 5.0, 2, "Schaffer F6" );
	print_benchmark<Rosenbrock>(opt, rosenbrock, init_4, 2, "Rosenbrock 2D" );
}

int main() {
	std::cout << std::setprecision(2) << std::scientific;

	// {
	// 	GradientDescent<Dynamic> optimizer(200, 1e-5, 1e-2);
	// 	test_optim(optimizer, "Gradient descent");
	// }

	// {
	// 	GradientDescent<Dynamic, WolfeLineSearch> optimizer(500, 1e-5, 1e-2);
	// 	test_optim(optimizer, "Gradient descent, Wolfe line search");
	// }

	{
		NelderMead<Dynamic> optimizer(500, 1e-5);
		test_optim(optimizer, "Nelder-Mead");
	}

	// {
	// 	ConjugateGradient<Dynamic, WolfeLineSearch> optimizer(500, 1e-5, 1e-2, true);
	// 	test_optim(optimizer, "Conjugate gradient, Wolfe line search");
	// }

	// {
	// 	LBFGS<Dynamic, WolfeLineSearch> optimizer(500, 1e-5, 1e-2, 30);
	// 	test_optim(optimizer, "LBFGS-30, Wolfe line search");
	// }

	return 0;
}
