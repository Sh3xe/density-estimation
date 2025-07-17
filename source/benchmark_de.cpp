#include "benchmark_de.hpp"
#include <fstream>
#include <sstream>

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > split_dataset( const Eigen::MatrixXd &dataset, size_t k, size_t i) {
	size_t slice_size = dataset.size() / k;
	size_t slice_index_begin = std::min(i*slice_size, dataset.size() - slice_size);

	Eigen::MatrixXd testing = dataset.block(slice_index_begin, 0, slice_size, dataset.cols());

	Eigen::MatrixXd training_1 = dataset.block(0, 0, slice_index_begin, dataset.cols());
	Eigen::MatrixXd training_2 = dataset.block(slice_index_begin + slice_size, 0, dataset.size() - slice_index_begin, dataset.cols());

	Eigen::MatrixXd training(training_1.rows()+training_2.rows(), training_1.cols());
	training << training_1, training_2;

	return std::make_pair(training, testing);
}

std::string file_name(
	const std::string &test_directory,
	const std::string &optimizer
) {
	std::stringstream ss;
	ss << std::setprecision(2);
	ss << "outputs/" << optimizer << "_" << test_directory << "_log_density.csv";
	return ss.str();
}

DETestScenario<Triangulation<2, 2>> load_test_2d(const std::string &dir_name) {
	// Geometry definition
	auto dir = path("./data") / path(dir_name);

	Triangulation<2, 2> space_discr(
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
	
	return DETestScenario<Triangulation<2, 2>>(dir_name, std::move(dataset), std::move(g_init), std::move(space_discr) );
}

DETestScenario<Triangulation<1, 1>> load_test_snp500(const std::string &dir_name) {
// Geometry definition
	auto dir = path("./data") / path(dir_name);

	Triangulation<1, 1> space_discr(-100.0, 100.0, 200);

	Eigen::VectorXd g_init = Eigen::VectorXd::Constant(200, 1, 1.0); // TODO: initial gaussian distribution
	
	// Data
	Eigen::MatrixXd dataset =
		read_csv<double>(dir / path("data_space.csv"))
		.as_matrix();
	
	return DETestScenario<Triangulation<1, 1>>(dir_name, std::move(dataset), std::move(g_init), std::move(space_discr) );
}

void print_benchmark_md( const std::vector<DEBenchmarkResult> &benchmark, std::ostream &os ) {
	os << "### " << benchmark[0].test_title << "\n";
	os << "Optimizer|Iterations|Loss|CV error|lambda|\n";
	os << "-|-|-|-|-\n";

	for(const auto &res: benchmark) {
		os << res.optimizer << "|" << res.iterations << "|" << res.loss << "|" << res.cv_error << "|" << res.lambda << "\n";
	}
}

void save_log_densities( const std::vector<DEBenchmarkResult> &benchmark ) {
	auto output_dir = path("outputs");
	for(const DEBenchmarkResult &res: benchmark) {
		std::string filename = output_dir / file_name(res.test_title, res.optimizer);
		write_csv(filename, res.log_density);
	}
}

void full_benchmark()
{
	// Load test scenario
	// auto snp500 = load_test_snp500("snp500");
	auto gaussian_square = load_test_2d("gaussian_square");
	auto infections_southampton = load_test_2d("infections_southampton");
	// auto horseshoe = load_test_scenario("horseshoe");
	// auto kent_sphere = load_test_scenario("kent_sphere");
	// auto gaussian_curved = load_test_scenario("gaussian_curved");
	// auto earthquake = load_test_scenario("earthquake");
	// auto accidents_bergamo = load_test_scenario("accidents_bergamo");

	// Optimizers
	auto lbfgs30 = std::make_pair(std::string{"lbfgs30"}, LBFGS<Eigen::Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, 30, ERR_TOL, STEP_SIZE});
	auto grad_descent = std::make_pair(std::string{"gd"}, GradientDescent<Eigen::Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE});
	// auto cg_pr = std::make_pair( std::string{""}, ConjugateGradient<Eigen::Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, true});
	// auto cg_fr = std::make_pair( std::string{""}, ConjugateGradient<Eigen::Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, false});
	// auto nelder_mead = std::make_pair( std::string{""}, NelderMead<Eigen::Dynamic> {ERR_TOL, STEP_SIZE});


	// Benchmarks
	std::fstream output_file {"outputs/benchmark.md", std::ios::out | std::ios::trunc};
	if(!output_file) {
		std::cerr << "Canont open \"outputs/benchmark.md\"" << std::endl;
		return;
	}

	{
		auto gaussian_square_res = benchmark_all_opt(gaussian_square, lbfgs30, grad_descent /*, cg_pr, cg_fr, nelder_mead*/);
		print_benchmark_md(gaussian_square_res, output_file);
		save_log_densities(gaussian_square_res);
	}

	{	
		auto infections_southampton_res = benchmark_all_opt(infections_southampton, lbfgs30, grad_descent /*, cg_pr, cg_fr, nelder_mead*/ );
		print_benchmark_md(infections_southampton_res, output_file);
	}
	// auto horseshoe_res = benchmark_all_opt(horseshoe, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(horseshoe_res, output_file);

	// auto kent_sphere_res = benchmark_all_opt(kent_sphere, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(kent_sphere_res, output_file);

	// auto gaussian_curved_res = benchmark_all_opt(gaussian_curved, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(gaussian_curved_res, output_file);

	// auto earthquake_res = benchmark_all_opt(earthquake, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(earthquake_res, output_file);

	// auto accidents_bergamo_res = benchmark_all_opt(accidents_bergamo, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(accidents_bergamo_res, output_file);
}