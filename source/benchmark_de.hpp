#pragma once

#include "utilities.hpp"
#include "pch.hpp"

namespace de {

constexpr double LINK_TOL = 1e-5;
constexpr int MAX_ITERATIONS = 200;
constexpr int CV_K = 5;
constexpr double STEP_SIZE = 1e-2;
constexpr double ERR_TOL = 1e-5;

/**
 * @brief Result of the benchmakr of a single (optimizer, test_case) pair
 * 
 */
struct DEBenchmarkResult {
	std::string test_title;
	std::string optimizer;

	size_t iterations;
	Eigen::VectorXd log_density;
	microsec duration;
	microsec cv_duration;

	double loss;
	double lambda;
	double cv_error;
};

template <typename SpaceTriangulation>
struct DETestScenario {

	static constexpr int embed_dim = SpaceTriangulation::embed_dim;

	std::string title;
	Eigen::MatrixXd dataset;
	Eigen::VectorXd g_init;
	SpaceTriangulation discretization;

	DETestScenario(const std::string &title, Eigen::MatrixXd &&dataset, Eigen::VectorXd &&g_init, SpaceTriangulation &&discretization):
		title( title),
		dataset( dataset),
		g_init( g_init),
		discretization( discretization) {
	}
};

/**
 * @brief Generate a pair of (Eigen::MatrixXd, test_set) of size (SIZE*(k-1/k), SIZE/k)
 * 
 * @param dataset original training set
 * @param k number of subdivisions of the training set
 * @param i which slice of the k will be used as test_set
 * @return std::pair< Eigen::MatrixXd, Eigen::MatrixXd > a pair (training, testing)
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > split_dataset( const Eigen::MatrixXd &dataset, size_t k, size_t i);

/**
 * @brief Generate a sequence of 10^x_i with x_i in [from, to]
 * 
 * @param from smaller than to, begin of range
 * @param to greater than from, end of range
 * @param incr increments
 * @return Eigen::VectorXd propositions for lambda
 */
Eigen::VectorXd gen_lambda_prop(double from, double to, double incr);

DETestScenario<fdapde::Triangulation<2, 2>> load_test_2d(const std::string &dir_name);

DETestScenario<fdapde::Triangulation<2, 3>> load_test_2_5d(const std::string &dir_name);

DETestScenario<fdapde::Triangulation<1, 2>> load_test_1_5_d(const std::string &dir_name);

DETestScenario<fdapde::Triangulation<1, 1>> load_test_snp500(const std::string &dir_name);

template <
	typename Optimizer,
	typename DEPDESolver,
	typename SpaceTriangulation,
	typename ...Hooks
>
double cv_error(
	Optimizer &optimizer,
	DEPDESolver &solver,
	DETestScenario<SpaceTriangulation> &scenario,
	double lambda,
	Hooks &&...hooks
) {
	double err_tot = 0.0;

	for(int i = 0; i < CV_K; ++i) {
		// Split the dataset & setup solver
		const auto [train_set, test_set] = split_dataset(scenario.dataset, CV_K, i);
		fdapde::GeoFrame geo_data(scenario.discretization);

		auto& sample = geo_data.template insert_scalar_layer<fdapde::POINT>("sample", train_set);
		fdapde::DEPDE<typename DEPDESolver::solver_t> model(geo_data, solver);
		model.set_llik_tolerance(LINK_TOL);

		// Solve
		model.fit(lambda, scenario.g_init, optimizer, std::forward<Hooks>(hooks)...);

		// Compute the error on the test set
		fdapde::FeSpace Vh(scenario.discretization, fdapde::P1<1>);
		fdapde::FeFunction log_dens_func(Vh, model.log_density());

		double err = 0.0;
		for(int i = 0; i < test_set.rows(); ++i) {
			err -= log_dens_func(test_set.row(i));
		}
		err_tot += err / static_cast<double>(test_set.rows());
	}

	return err_tot / static_cast<double>(CV_K);
}

template <
	typename Optimizer,
	typename SpaceTriangulation,
	typename ...Hooks
>
DEBenchmarkResult benchmark_one(
	Optimizer &optimizer,
	DETestScenario<SpaceTriangulation> &scenario,
	const Eigen::VectorXd &lambda_prop,
	const std::string &optimizer_title,
	Hooks &&...hooks
) {
	DEBenchmarkResult res;
	res.test_title = scenario.title;
	res.optimizer = optimizer_title;

	// Creating the solver
	fdapde::FeSpace Vh(scenario.discretization, fdapde::P1<1>);
	fdapde::TrialFunction f(Vh);
	fdapde::TestFunction v(Vh);
	auto a = fdapde::integral(scenario.discretization)(dot(grad(f), grad(v)));
	fdapde::ZeroField<DETestScenario<SpaceTriangulation>::embed_dim> u;
	auto F = fdapde::integral(scenario.discretization)(u * v);
	auto solver = fdapde::fe_de_elliptic(a, F);
	
	// Find lambda that minimized the CV error
	res.cv_error = std::numeric_limits<double>::max();
	res.lambda = lambda_prop[0];
	auto cv_begin_time = std::chrono::high_resolution_clock::now();
	for(int i = 0; i < lambda_prop.rows(); ++i) {
		std::cout << "Computing CV err (" << i+1 << "/" << lambda_prop.rows() <<") for " << scenario.title << " " << optimizer_title << std::endl;
		double err = cv_error(optimizer, solver, scenario, lambda_prop[i], std::forward<Hooks>(hooks)...);
		if(err < res.cv_error) {
			res.cv_error = err;
			res.lambda = lambda_prop[i];
		}
	}
	res.cv_duration = std::chrono::high_resolution_clock::now() - cv_begin_time;
	
	// Now that we know the à-priori best lambda value, let's do the full benchmark
	fdapde::GeoFrame geo_data(scenario.discretization);
	auto& sample = geo_data.template insert_scalar_layer<fdapde::POINT>("sample", scenario.dataset);
	fdapde::DEPDE model(geo_data, solver);
	model.set_llik_tolerance(LINK_TOL);
	auto begin_time = std::chrono::high_resolution_clock::now();
	model.fit(res.lambda, scenario.g_init, optimizer, std::forward<Hooks>(hooks)...);
	res.duration = std::chrono::high_resolution_clock::now() - begin_time;
	res.log_density = model.log_density();
	res.iterations = optimizer.n_iter();
	res.loss = optimizer.value();

	return res;
}

template <typename DETestScenarioType>
std::vector<DEBenchmarkResult> benchmark_all_opt(DETestScenarioType &scenario, const Eigen::VectorXd &lambda_prop) {
	// Benchmark one optimizer
	std::vector<DEBenchmarkResult> results;
	
	std::cout << "Starting benchmark for " << scenario.title << std::endl;

	{
		auto lbfgs30 = fdapde::LBFGS<fdapde::Dynamic> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, 30};
		results.push_back(benchmark_one(
			lbfgs30, scenario, lambda_prop, "lbfgs30",
			fdapde::WolfeLineSearch()
		));
		std::cout << scenario.title << ": lbfgs30 done" << std::endl;
	}

	{
		auto grad_descent = fdapde::GradientDescent<fdapde::Dynamic> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE};
		results.push_back(benchmark_one(
			grad_descent, scenario, lambda_prop,"grad_descent",
			fdapde::WolfeLineSearch()
		));
		std::cout << scenario.title << ": grad_descent done" << std::endl;
	}

	{
		auto cg_pr = fdapde::PolakRibiereCG<fdapde::Dynamic> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE};
		results.push_back(benchmark_one(
			cg_pr, scenario, lambda_prop, "cg_pr",
			fdapde::WolfeLineSearch()
		));
		std::cout << scenario.title << ": cg_pr done" << std::endl;
	}

	{
		auto cg_fr = fdapde::FletcherReevesCG<fdapde::Dynamic> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE};
		results.push_back(benchmark_one(
			cg_fr, scenario, lambda_prop, "cg_fr",
			fdapde::WolfeLineSearch()
		));
		std::cout << scenario.title << ": cg_fr done" << std::endl;
	}

	return results;
}

void print_benchmark_md( const std::vector<DEBenchmarkResult> &benchmark, std::ostream &os );

void save_log_densities( const std::vector<DEBenchmarkResult> &benchmark );

void full_benchmark(bool output_csv = false);

} // namespace de