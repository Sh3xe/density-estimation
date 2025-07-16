#include <fdaPDE/fdapde.h>
#include <filesystem>
#include <sstream>
#include <string>
#include <chrono>
#include <cstdint>
#include <utility>
#include <ostream>

using namespace fdapde;
using path = std::filesystem::path;
using microsec = std::chrono::duration<double, std::micro>;
using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
using vector_t = Eigen::Matrix<double, Dynamic, 1>;

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
	vector_t log_density;
	microsec duration;
	mecrosec avg_cv_duration;

	double loss;
	double lambda;
	double cv_error;
};

template <typename SpaceTriangulation>
struct DETestScenario {
	DETestScenario() = default;

	std::string &title;
	matrix_t dataset;
	vector_t g_init;
	SpaceTriangulation discretization;
}

/**
 * @brief Generate a pair of (matrix_t, test_set) of size (SIZE*(k-1/k), SIZE/k)
 * 
 * @param dataset original training set
 * @param k number of subdivisions of the training set
 * @param i which slice of the k will be used as test_set
 * @return std::pair< matrix_t, matrix_t > a pair (training, testing)
 */
std::pair< matrix_t, matrix_t > split_dataset( const matrix_t &dataset, size_t k, size_t i) {
	size_t slice_size = dataset.size() / k;
	size_t slice_index_begin = std::min(i*slice_size, dataset.size() - slice_size);

	matrix_t testing = dataset.block(slice_index_begin, 0, slice_size, dataset.cols());

	matrix_t training_1 = dataset.block(0, 0, slice_index_begin, dataset.cols());
	matrix_t training_2 = dataset.block(slice_index_begin + slice_size, 0, dataset.size() - slice_index_begin, dataset.cols());

	matrix_t training(training_1.rows()+training_2.rows(), training_1.cols());
	training << training_1, training_2;

	return std::make_pair(training, testing);
}

/**
 * @brief Construct a file path for the output log density given some parameters
 * 
 * @param test_directory 
 * @param optimizer 
 * @return std::string 
 */
std::string file_name(
	const std::string &test_directory,
	const std::string &optimizer
) {
	std::stringstream ss;
	ss << std::setprecision(2);
	ss << "outputs/" << optimizer << "_" << test_directory << "_log_density.csv";
	return ss.str();
}

template <typename SpaceTriangulation>
DETestScenario<SpaceTriangulation> load_test_scenario(const std::string &dir_name) {
	// Geometry definition
	auto dir = path("./data") / path(dir_name);

	Triangulation<2, 2> space_discr(
		(dir / path("mesh_vertices.csv")).string(),
		(dir / path("mesh_elements.csv")).string(),
		(dir / path("mesh_boundary.csv")).string(),
		true, true
	);

	vector_t g_init =
	read_csv<double>(dir / path("f_init.csv"))
		.as_matrix()
		.array()
		.log();
	
	// Data
	matrix_t dataset =
		read_csv<double>(dir / path("sample.csv"))
		.as_matrix();
	
	return DETestScenario( dir_name, dataset, g_init, space_discr );
}

template <typename SpaceTriangulation>
DETestScenario<SpaceTriangulation> load_test_snp500(const std::string &dir_name) {
	// Geometry definition
	auto dir = path("./data") / path(dir_name);

	Triangulation<1, 1> space_discr(-100.0, 100.0, 200);

	vector_t g_init = Eigen::Constant(200, 1, 1.0); // TODO: initial gaussian distribution
	
	// Data
	matrix_t dataset =
		read_csv<double>(dir / path("data_space.csv"))
		.as_matrix();
	
	return DETestScenario( dir_name, dataset, g_init, space_discr );
}

template <
	typename Optimizer,
	typename DEPDESolver,
	typename SpaceTriangulation
>
DEBenchmarkResult benchmark_one(
	Optimizer &optimizer,
	DEPDESolver &solver,
	const DETestScenario<SpaceTriangulation> &scenario
	// vector_t lambda_prop,
	const std::string &optimizer_title
) {
	DEBenchmarkResult res;
	
	// calculate the cross_validation error
	res.cv_error = 0.0;
	res.avg_cv_duration = microsec(0);
	for(int i = 0; i < CV_K; ++i) {
		// Split the dataset & setup solver
		const auto [train_set, test_set] = split_dataset(scenario.dataset, CV_K, i);
		GeoFrame geo_data(scenario.discretization);
		auto& l1 = geo_data.insert_scalar_layer<POINT>("l1", train_set);
		DEPDE m(geo_data, solver);
		m.set_llik_tolerance(LINK_TOL);

		// Solve
		auto begin_time = std::chrono::high_resolution_clock::now();
		solver.fit(0.002, scenario.g_init, optimizer);
		auto end_time = std::chrono::high_resolution_clock::now();
		res.avg_cv_duration += end_time - begin_time;

		// Compute the error on the test set
		res.cv_error += optimizer.value();
		// TODO: How to compute the log likelyhood given the space discretization and the log-density ?
	}

	res.lambda = 0.002; // TODO: find the best lambda value
	res.cv_error /= static_cast<double>(k);
	res.avg_cv_duration /= static_cast<double>(k);
	
	res.test_title = scenario.title;
	res.optimizer = optimizer_title;
	
	// Now that we know the à-priori best lambda value, let's do the full benchmark
	GeoFrame geo_data(scenario.discretization);
	auto& l1 = geo_data.insert_scalar_layer<POINT>("l1", dataset);
	DEPDE m(geo_data, solver);
	m.set_llik_tolerance(LINK_TOL);
	auto begin_time = std::chrono::high_resolution_clock::now();
	solver.fit(res.lambda, scenario.g_init, optimizer);
	auto end_time = std::chrono::high_resolution_clock::now();
	res.duration = end_time - begin_time;
	res.log_density = solver.log_density();
	res.iterations = optimizer.n_iter();
	res.loss = optimizer.value();

	return res;
}

template <typename DETestScenario, typename Solver, typename ...NamedOptimizers>
std::vector<DEBenchmarkResult> benchmark_all_opt(const DETestScenario &test_scenario, const Solver &solver, NamedOptimizers &opt_head, NamedOptimizers &...opt_tail) {
	// Benchmark one optimizer
	auto &[opt_name, optimizer] = opt_head;
	DEBenchmarkResult benchmark_res = benchmark_one( optimizer, solver, test_scenario, opt_name );

	// Compile-time-recursively benchmark the rest
	auto result_rest = benchmark_all_opt(test_scenario, solver, ...opt_tail);

	// Combine the results into a single list
	std::vector<DEBenchmarkResult> results;
	results.insert(results.end(), result_rest.begin(), result_rest.end());
	return results;
}

/**
 * @brief Print to 
 * 
 * @param benchmark 
 * @param os 
 */
void print_benchmark_md( const std::vector<DEBenchmarkResult> &benchmark, std::ostream &os );

/**
 * @brief Perform the full benchmark, i.e. benchmakr CV error, for each pair (optimizer, test_case)
 * 
 */
void full_benchmark();