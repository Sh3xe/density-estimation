#include <fdaPDE/optimization.h>
#include <iomanip>
#include <chrono>
#include <Eigen/Dense>

using namespace fdapde;
using vector_t = Eigen::Matrix<double, Dynamic, 1>;
using microsec = std::chrono::duration<double, std::micro>;

struct SchafferF6 {
	SchafferF6() = default;

	static vector_t argmin(int dim) {
		return vector_t::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const vector_t& x) const {
		double r = x[0]*x[0] + x[1]*x[1];
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
		return 418.9829 + dimension * s;
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

struct SquareFunction {
	SquareFunction() = default;

	static vector_t argmin(int dim) {
		return vector_t::Zero(dim, 1);
	}

	static double min() { return 0.0; }

	double operator()(const vector_t& x) const {
		double res = 0.0;

		for(int i = 0; i < x.rows(); ++i)
			res += x[i]*x[i];

		return res;
	}
};

template <typename OptimizerType>
void print_performances(
	OptimizerType &opt,
	const vector_t &argmin,
	double min,
	double time_microsec,
	const std::string &test_case
) {
	// print a markdown table, to be exported to pdf for better visibility
	double x_err = (argmin - opt.optimum()).norm();
	double f_err = (min - opt.value()); f_err *= f_err;
	std::cout << 
		"|" << test_case << 
		"|" << opt.n_iter() << 
		"|" << x_err <<
		"|" << f_err <<
		"|" << min << 
		"|" << opt.value() << 
		"|" << time_microsec << 
		"|\n";
}

template<
	typename OptimizerFunc,
	typename Optimizer,
	typename FieldFunc
>
void test_optimizer_to_md(
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
		duration / 100.0,
		name
	);
}

template <typename OptimizerType>
void test_optim(OptimizerType &opt, const std::string &title) {

	// Optimization functions and their init point
	ScalarField<Dynamic, SchafferF6>     schaffer_f6(2);
	ScalarField<Dynamic, Schwefel>       schwefel_2d(2);
	ScalarField<Dynamic, Schwefel>       schwefel_10d(10);
	ScalarField<Dynamic, SquareFunction> square_2d(2);
	ScalarField<Dynamic, SquareFunction> square_30d(30);
	ScalarField<Dynamic, Rastrigin>      rastrigin_30d(30);
	ScalarField<Dynamic, Rastrigin>      rastrigin_2d(2);
	
	// Initial points
	vector_t init_1 = vector_t::Constant(2, 1, 1.0);
	vector_t init_2 = vector_t::Constant(10, 1, 1.0);
	vector_t init_3 = vector_t::Constant(30, 1, 1.0);

	// Markdown header
	std::cout << "### `" << title << "`\n";
	std::cout << "|Method|n_iter| x-x* l2 err | f-f* l2 err | real min | approx min | time (microsec) |" << std::endl;
	std::cout << "|-|-|-|-|-|-|-|" << std::endl;

	test_optimizer_to_md<SquareFunction>(opt, square_2d, init_1, 2, "Square 2D" );
	test_optimizer_to_md<SquareFunction>(opt, square_30d, init_3, 30, "Square 30D" );

	test_optimizer_to_md<Schwefel>(opt, schwefel_2d, init_1 * 20.0, 2, "Schwefel 2D" );
	test_optimizer_to_md<Schwefel>(opt, schwefel_10d, init_2 * 20.0, 10, "Schwefel 10D" );

	test_optimizer_to_md<Rastrigin>(opt, rastrigin_2d, init_1 * 3.0, 2, "Rastrigin 2D" );
	test_optimizer_to_md<Rastrigin>(opt, rastrigin_30d, init_3 * 3.0, 30, "Rastrigin 30D" );

	test_optimizer_to_md<SchafferF6>(opt, schaffer_f6, init_1, 2, "Schaffer F6" );
}

int main() {
	std::cout << std::setprecision(2);

	{
		GradientDescent<Dynamic, WolfeLineSearch> optimizer(500, 1e-5, 1e-2);
		test_optim(optimizer, "GradientDescent<Dynamic, WolfeLineSearch>");
	}

	return 0;
}
