#include "benchmark_de.hpp"
#include "utilities.hpp"
#include "cmake_defines.hpp"

#include <fdaPDE/fdapde.h>

#include <filesystem>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>
#include <utility>
#include <thread>
#include <mutex>

using namespace fdapde;
namespace fs = std::filesystem;
using namespace Eigen;

constexpr double LINK_TOL = 1e-5;
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
 * @brief Generate a pair of (train_set, test_set) of size (SIZE*(k-1/k), SIZE/k)
 * 
 * @param dataset original training set
 * @param k number of subdivisions of the training set
 * @param i which slice of the k will be used as test_set
 * @return std::pair< Eigen::MatrixXd, Eigen::MatrixXd > a pair (training, testing)
 */
std::pair< MatrixXd, MatrixXd > split_dataset( const MatrixXd &dataset, size_t k, size_t i) {
	size_t slice_size = dataset.rows() / k;
	size_t slice_index_begin = std::min(i*slice_size, dataset.rows() - slice_size);

	MatrixXd testing = dataset.block(slice_index_begin, 0, slice_size, dataset.cols());

	MatrixXd training_2 = dataset.block(slice_index_begin + slice_size, 0, dataset.rows() - slice_index_begin - slice_size, dataset.cols());
	MatrixXd training_1 = dataset.block(0, 0, slice_index_begin, dataset.cols());

	MatrixXd training(training_1.rows()+training_2.rows(), training_1.cols());
	training << training_1, training_2;

	return std::make_pair(training, testing);
}

/**
 * @brief Generate a sequence of 10^x_i with x_i in [from, to]
 * 
 * @param from smaller than to, begin of range
 * @param to greater than from, end of range
 * @param incr increments
 * @return Eigen::VectorXd propositions for lambda
 */
VectorXd gen_lambda_prop(double from, double to, double incr) {
	std::vector<double> props;
	for(double exponent = from; exponent <= to + 1e-4; exponent += incr) {
		props.push_back(std::pow(10.0, exponent));
	}

	VectorXd lambda_prop(props.size());
	for(size_t i = 0; i < props.size(); ++i)
		lambda_prop[i] = props[i];

	return lambda_prop;
}

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
	const std::string &opt_name,
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

		write_csv("cv_csv/cpp#" + std::to_string(i) + "#" + opt_name + "#" + std::to_string(lambda) + "#" + scenario.title + ".csv", model.log_density());

		double total_dens_eval = 0.0;
		for(int i = 0; i < test_set.rows(); ++i) {
			total_dens_eval += std::exp( log_dens_func(test_set.row(i)) );
		}

		err_tot += fdapde::integral(scenario.discretization, fdapde::QS2DP4)(exp(2.0*log_dens_func));
		err_tot	-= 2.0 * (total_dens_eval / (double)test_set.rows());
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
		double err = cv_error(optimizer, solver, scenario, lambda_prop[i], optimizer_title, std::forward<Hooks>(hooks)...);
		std::cout << "CV err for " << lambda_prop[i]<< " " << optimizer_title << " = " << err << ", " << scenario.title << std::endl;
		if(err < res.cv_error && !std::isnan(err)) {
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
		auto bfgs = fdapde::BFGS<fdapde::Dynamic> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE};
		results.push_back(benchmark_one(
			bfgs, scenario, lambda_prop,"BFGS",
			fdapde::WolfeLineSearch()
		));
		std::cout << scenario.title << ": BFGS done" << std::endl;
	}

	{
		auto lbfgs30 = fdapde::LBFGS<fdapde::Dynamic> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, 10};
		results.push_back(benchmark_one(
			lbfgs30, scenario, lambda_prop, "LBFGS10",
			fdapde::WolfeLineSearch()
		));
		std::cout << scenario.title << ": LBFGS10 done" << std::endl;
	}

	{
		auto cg_pr = fdapde::PolakRibiereCG<fdapde::Dynamic> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, true};
		results.push_back(benchmark_one(
			cg_pr, scenario, lambda_prop, "cg_pr",
			fdapde::WolfeLineSearch()
		));
		std::cout << scenario.title << ": cg_pr done" << std::endl;
	}

	return results;
}

void save_log_densities( const std::vector<DEBenchmarkResult> &benchmark ) {
	auto output_dir = path(ROOT_DIR) / path("outputs");
	for(const DEBenchmarkResult &res: benchmark) {
		path filename = output_dir / path("cpp_" + utils::file_name(res.test_title, res.optimizer));
		std::cout << "filename" << filename.string() << std::endl;
		write_csv(filename.string(), res.log_density);
	}
}

void print_benchmark_md( const std::vector<DEBenchmarkResult> &benchmark, std::ostream &os ) {
	static std::mutex lock;

	lock.lock();
	os << "### " << benchmark[0].test_title << "\n";
	os << "Optimizer|Iterations|Loss|CV error|lambda|cv_duration|duration|\n";
	os << "-|-|-|-|-|-|-\n";

	for(const auto &res: benchmark) {
		os
			<< res.optimizer << "|"
			<< res.iterations << "|"
			<< res.loss << "|"
			<< res.cv_error << "|"
			<< res.lambda << "|"
			<< utils::duration_to_str(res.cv_duration) << "|"
			<< utils::duration_to_str(res.duration) << "|"
			<< "\n";
	}

	os << "\n";

	lock.unlock();
}

// Benchmark

DETestScenario<Triangulation<2, 2>> load_test_2d(const std::string &dir_name) {
	std::cout << "Loading 2D " << dir_name << std::endl;
	// Geometry definition
	auto dir = path(ROOT_DIR) / path("data") / path(dir_name);

	Triangulation<2, 2> space_discr(
		(dir / path("mesh_vertices.csv")).string(),
		(dir / path("mesh_elements.csv")).string(),
		(dir / path("mesh_boundary.csv")).string(),
		true, true
	);

	Eigen::MatrixXd dataset =
		read_csv<double>(dir / path("sample.csv"))
		.as_matrix();

	Eigen::VectorXd g_init =
		read_csv<double>(dir / path("f_init.csv"))
		.as_matrix()
		.array()
		.log();

	std::cout << "Done loading 2D " << dir_name << std::endl;
	
	return DETestScenario<Triangulation<2, 2>>(dir_name, std::move(dataset), std::move(g_init), std::move(space_discr) );
}

DETestScenario<Triangulation<2, 3>> load_test_2_5d(const std::string &dir_name) {
	std::cout << "Loading 2.5D " << dir_name << std::endl;
	// Geometry definition
	auto dir = path(ROOT_DIR) / path("data") / path(dir_name);

	Triangulation<2, 3> space_discr(
		(dir / path("mesh_vertices.csv")).string(),
		(dir / path("mesh_elements.csv")).string(),
		(dir / path("mesh_boundary.csv")).string(),
		true, true
	);

	Eigen::VectorXd g_init =
		read_csv<double>(dir / path("f_init.csv"))
		.as_matrix()
		.array()
		.log();

	// Data
	Eigen::MatrixXd dataset =
		read_csv<double>(dir / path("sample.csv"))
		.as_matrix();

	std::cout << "Done loading 2.5D " << dir_name << std::endl;
	
	return DETestScenario<Triangulation<2, 3>>(dir_name, std::move(dataset), std::move(g_init), std::move(space_discr) );
}

DETestScenario<Triangulation<1, 2>> load_test_1_5_d(const std::string &dir_name) {
	std::cout << "Loading 1.5D " << dir_name << std::endl;
	// Geometry definition
	auto dir = path(ROOT_DIR) / path("data") / path(dir_name);

	Triangulation<1, 2> space_discr(
		(dir / path("mesh_vertices.csv")).string(),
		(dir / path("mesh_edges.csv")).string(),
		(dir / path("mesh_boundary.csv")).string(),
		true, true
	);
	
	Eigen::VectorXd g_init =
		read_csv<double>(dir / path("f_init.csv"))
		.as_matrix()
		.array()
		.log();

	// Data
	Eigen::MatrixXd dataset =
		read_csv<double>(dir / path("sample.csv"))
		.as_matrix();

	std::cout << "Done loading 1.5D " << dir_name << std::endl;
	
	return DETestScenario<Triangulation<1, 2>>(dir_name, std::move(dataset), std::move(g_init), std::move(space_discr) );
}


template <
	typename SpaceTriangulation
>
double l2_error(
	DETestScenario<SpaceTriangulation> &scenario,
	const Eigen::VectorXd &true_density_vec, 
	const Eigen::VectorXd &estimated_density_vec
) {

	fdapde::FeSpace Vh(scenario.discretization, fdapde::P1<1>);
	fdapde::FeFunction diff_func(Vh, (true_density_vec - estimated_density_vec));

	return fdapde::integral(scenario.discretization, fdapde::QS2DP4)( diff_func * diff_func );
}

void print_l2_errors(const std::string &test_path) {
	// fetch real density in the data directory
	auto dir = path(ROOT_DIR) / path("data") / path(test_path);
	auto scenario = load_test_2d(test_path);

	auto true_density_vec = read_csv<double>(dir / path("true_density.csv"))
		.as_matrix()
		.array();
		std::cout << true_density_vec.rows() << " " << true_density_vec.cols() << std::endl;

	for( auto &file: fs::directory_iterator(path(ROOT_DIR) / path("outputs")) ){
		// Filter csv files
		auto filename = file.path().filename().string();
		std::cout << "1" << filename << std::endl;
		if( file.path().extension() != ".csv") continue;	

		std::cout << "2" << std::endl;

		// Filter for the right test case
		if( filename.find(test_path) == std::string::npos ) continue;
	
		std::cout << "3 " << filename << std::endl;
		std::cout << "3 " << file.path().string() << std::endl;

		// Compute L2 error
		auto estimated_density_vec = read_csv<double>(file.path().string())
			.as_matrix()
			.array()
			.exp();

		std::cout << "4" << std::endl;
		
		std::cout << estimated_density_vec.rows() << " " << estimated_density_vec.cols() << std::endl;

		double err = l2_error(scenario, true_density_vec, estimated_density_vec);
		std::cout << filename << ": " << err << std::endl;
	}
}


void de_full_benchmark(bool output_csv)
{
	// Output markdown file
	std::fstream output_file {path(ROOT_DIR) / path("outputs/cpp_de_bmrk.md"), std::ios::out | std::ios::trunc};
	if(!output_file) {
		std::cerr << "Cannot open \"outputs/benchmark.md\"" << std::endl;
		return;
	}

	// Creating one thread per benchmark
	std::vector<std::thread> workers;
	auto lambda_proposal = gen_lambda_prop(-2.0, 2.0, 1.0);

	workers.emplace_back( [&](){
		auto gaussian_square = load_test_2d("gaussian_square");
		auto gaussian_square_res = benchmark_all_opt(gaussian_square, lambda_proposal);
		print_benchmark_md(gaussian_square_res, output_file);
		if(output_csv) save_log_densities(gaussian_square_res);
	});

	workers.emplace_back([&](){	
		auto infections_southampton = load_test_2d("infections_southampton");
		auto infections_southampton_res = benchmark_all_opt(infections_southampton, lambda_proposal);
		print_benchmark_md(infections_southampton_res, output_file);
		if(output_csv) save_log_densities(infections_southampton_res);
	});

	workers.emplace_back([&](){
		auto horseshoe = load_test_2d("horseshoe");
		auto horseshoe_res = benchmark_all_opt(horseshoe, lambda_proposal);
		print_benchmark_md(horseshoe_res, output_file);
		if(output_csv) save_log_densities(horseshoe_res);
	});

	// workers.emplace_back([&](){
	// 	auto kent_sphere = load_test_2_5d("kent_sphere");
	// 	auto kent_sphere_res = benchmark_all_opt(kent_sphere, gen_lambda_prop(-5.0, -1.0, 1.0));
	// 	print_benchmark_md(kent_sphere_res, output_file);
	// 	if(output_csv)
	// 	save_log_densities(kent_sphere_res);
	// });

	// workers.emplace_back([&](){
	// 	auto curved = load_test_2_5d("curved");
	// 	auto curved_res = benchmark_all_opt(curved, gen_lambda_prop(-4.0, 0.0, 1.0));
	// 	print_benchmark_md(curved_res, output_file);
	// 	if(output_csv)
	// 	save_log_densities(curved_res);
	// });

	// workers.emplace_back([&](){
	// 	auto earthquake = load_test_2_5d("earthquake");
	// 	auto earthquake_res = benchmark_all_opt(earthquake, Eigen::VectorXd::Constant(1, 1, 0.1));
	// 	print_benchmark_md(earthquake_res, output_file);
	// 	if(output_csv)
	// 	save_log_densities(earthquake_res);
	// });

	// workers.emplace_back([&](){
		// auto accidents_bergamo_res = benchmark_all_opt(accidents_bergamo);
		// print_benchmark_md(accidents_bergamo_res, output_file);
		// if(output_csv)
		// 	save_log_densities(accidents_bergamo_res)
	// });

	// workers.emplace_back([&](){
	// 	auto snp500 = load_test_snp500("snp500");
	// 	auto snp500_res = benchmark_all_opt(snp500, gen_lambda_prop(-2.0, 2.0, 1.0));
	// 	print_benchmark_md(snp500_res, output_file);
	// 	if(output_csv) save_log_densities(snp500_res);
	// });

	for(auto &t: workers) {
		t.join();
	}
}