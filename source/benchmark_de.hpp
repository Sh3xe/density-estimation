#include <fdaPDE/fdapde.h>
#include <filesystem>
#include <string>
#include <chrono>
#include <utility>

using namespace fdapde;
using path = std::filesystem::path;
using microsec = std::chrono::duration<double, std::micro>;

template <typename OptimizerType>
using NamedOptimizer = std::pair<std::string, OptimizerType>;

constexpr double LINK_TOL = 1e-8;
constexpr int MAX_ITERATIONS = 500;
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
	microsec avg_cv_duration;

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
 * @brief Construct a file path for the output log density given some parameters
 * 
 * @param test_directory 
 * @param optimizer 
 * @return std::string 
 */
std::string file_name( const std::string &test_directory, const std::string &optimizer );

DETestScenario<Triangulation<2, 2>> load_test_2d(const std::string &dir_name);

DETestScenario<Triangulation<2, 3>> load_test_2_5d(const std::string &dir_name);

DETestScenario<Triangulation<1, 1>> load_test_snp500(const std::string &dir_name);

template <
	typename Optimizer,
	typename DEPDESolver,
	typename SpaceTriangulation
>
DEBenchmarkResult benchmark_one(
	Optimizer &optimizer,
	DEPDESolver &solver,
	DETestScenario<SpaceTriangulation> &scenario,
	// Eigen::VectorXd lambda_prop,
	const std::string &optimizer_title
) {
	DEBenchmarkResult res;
	
	// calculate the cross_validation error
	// res.cv_error = 0.0;
	// res.avg_cv_duration = microsec(0);
	// for(int i = 0; i < CV_K; ++i) {
	// 	// Split the dataset & setup solver
	// 	const auto [train_set, test_set] = split_dataset(scenario.dataset, CV_K, i);
	// 	GeoFrame geo_data(scenario.discretization);
	// 	auto& l1 = geo_data.insert_scalar_layer<POINT>("l1", train_set);
	// 	DEPDE m(geo_data, solver);
	// 	m.set_llik_tolerance(LINK_TOL);

	// 	// Solve
	// 	auto begin_time = std::chrono::high_resolution_clock::now();
	// 	solver.fit(0.002, scenario.g_init, optimizer);
	// 	auto end_time = std::chrono::high_resolution_clock::now();
	// 	res.avg_cv_duration += end_time - begin_time;

	// 	// Compute the error on the test set
	// 	double error = 0.0;
	// 	// TODO: How to compute the log likelyhood given the space discretization and the log-density ?
	// }

	res.lambda = 0.02; // TODO: find the best lambda value
	// res.cv_error /= static_cast<double>(CV_K);
	// res.avg_cv_duration /= static_cast<double>(CV_K);
	
	res.test_title = scenario.title;
	res.optimizer = optimizer_title;
	
	// Now that we know the à-priori best lambda value, let's do the full benchmark

	GeoFrame geo_data(scenario.discretization);
	auto& sample = geo_data.template insert_scalar_layer<POINT>("sample", scenario.dataset);
	DEPDE<typename DEPDESolver::solver_t> model(geo_data, solver);
	model.set_llik_tolerance(LINK_TOL);
	auto begin_time = std::chrono::high_resolution_clock::now();
	model.fit(res.lambda, scenario.g_init, optimizer);
	auto end_time = std::chrono::high_resolution_clock::now();
	res.duration = end_time - begin_time;
	res.log_density = model.log_density();
	res.iterations = optimizer.n_iter();
	res.loss = optimizer.value();

	return res;
}

template <typename DETestScenarioType>
std::vector<DEBenchmarkResult> benchmark_all_opt(DETestScenarioType &scenario) {
	// Create solver
	FeSpace Vh(scenario.discretization, P1<1>);
	TrialFunction f(Vh);
	TestFunction  v(Vh);
	auto a = integral(scenario.discretization)(dot(grad(f), grad(v)));
	ZeroField<DETestScenarioType::embed_dim> u;
	auto F = integral(scenario.discretization)(u * v);
	auto solver = fe_de_elliptic(a, F);

	// Benchmark one optimizer
	std::vector<DEBenchmarkResult> results;

	auto lbfgs30 = LBFGS<Eigen::Dynamic, WolfeLineSearch> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, 30};
	auto grad_descent = GradientDescent<Eigen::Dynamic, WolfeLineSearch> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE};
	// auto cg_pr = ConjugateGradient<Eigen::Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, true});
	// auto cg_fr = ConjugateGradient<Eigen::Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, false});
	// auto nelder_mead = NelderMead<Eigen::Dynamic> {ERR_TOL, STEP_SIZE});

	std::cout << "Starting benchmark for " << scenario.title << std::endl;
	results.push_back(benchmark_one(lbfgs30, solver, scenario, "lbfgs30"));
	std::cout << scenario.title << ": lbfgs30 done" << std::endl;
	results.push_back(benchmark_one(grad_descent, solver, scenario, "grad_descent"));
	std::cout << scenario.title << ": grad_descent done" << std::endl;
	// results.push_back(benchmark_one(cg_pr, solver, scenario, "cg_pr"));
	// results.push_back(benchmark_one(cg_fr, solver, scenario, "cg_fr"));
	// results.push_back(benchmark_one(nelder_mead, solver, scenario, "nelder_mead"));

	return results;
}

void print_benchmark_md( const std::vector<DEBenchmarkResult> &benchmark, std::ostream &os );

void save_log_densities( const std::vector<DEBenchmarkResult> &benchmark );

void full_benchmark();