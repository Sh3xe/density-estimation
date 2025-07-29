#include "benchmark_de.hpp"

using namespace Eigen;
using namespace fdapde;

namespace de {

// Utilites

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

VectorXd gen_lambda_prop(double from, double to, double incr) {
	std::vector<double> props;
	for(double exponent = from; exponent < to; exponent += incr) {
		props.push_back(std::pow(10.0, exponent));
	}

	VectorXd lambda_prop(props.size());
	for(size_t i = 0; i < props.size(); ++i)
		lambda_prop[i] = props[i];

	return lambda_prop;
}

// Input-Output

void save_log_densities( const std::vector<DEBenchmarkResult> &benchmark ) {
	auto output_dir = path("outputs");
	for(const DEBenchmarkResult &res: benchmark) {
		std::string filename = output_dir / utils::file_name(res.test_title, res.optimizer);
		std::cout << "filename" << filename << std::endl;
		write_csv(filename, res.log_density);
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

	ScalarField<2, decltype([](const Eigen::VectorXd& p) { return 1.0; })> constant_unit;
	double area = integral(space_discr, QS2DP4)(constant_unit);
	Eigen::VectorXd g_init = Eigen::VectorXd::Constant(space_discr.nodes().rows(), 1, std::log(1.0/area));
	
	// Data
	Eigen::MatrixXd dataset =
		read_csv<double>(dir / path("sample.csv"))
		.as_matrix();

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

	// ScalarField<2, decltype([](const Eigen::VectorXd& p) { return 1.0; })> constant_unit;
	// double area = integral(space_discr, QS2DP4)(constant_unit);
	// Eigen::VectorXd g_init = Eigen::VectorXd::Constant(space_discr.nodes().rows(), 1, std::log(1.0/area));
	
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
	
	// Compute the integral over the whole domain for f_init
	double total_length = 0.0;
	for(int i = 0; i < space_discr.cells().rows(); ++i) {
		const VectorXi &cell = space_discr.cells().row(i);
		const auto &p1 = space_discr.nodes().row(cell[0]);
		const auto &p2 = space_discr.nodes().row(cell[1]);
		total_length += (p2 - p1).norm();
	}
	Eigen::VectorXd g_init = Eigen::VectorXd::Constant(space_discr.nodes().rows(), 1, std::log(1.0/total_length));

	// Data
	Eigen::MatrixXd dataset =
		read_csv<double>(dir / path("sample.csv"))
		.as_matrix();

	std::cout << "Done loading 1.5D " << dir_name << std::endl;
	
	return DETestScenario<Triangulation<1, 2>>(dir_name, std::move(dataset), std::move(g_init), std::move(space_discr) );
}

DETestScenario<Triangulation<1, 1>> load_test_snp500(const std::string &dir_name) {
	auto dir = path(ROOT_DIR) / path("data") / path(dir_name);
	
	// Data loading
	Eigen::VectorXd dataset = read_csv<double>( dir / path("data_space.csv"), true, true).as_matrix();
	double min = dataset.minCoeff();
	double max = dataset.maxCoeff();	
	constexpr int subdivisions = 300;
	
	// Geometry definition
	Triangulation<1, 1> space_discr(min - 1.0, max + 1.0, subdivisions);

	Eigen::VectorXd g_init = Eigen::VectorXd::Constant(subdivisions, 1, 1.0 / static_cast<double>(subdivisions));

	return DETestScenario<Triangulation<1, 1>>(dir_name, std::move(dataset), std::move(g_init), std::move(space_discr) );
}

void full_benchmark(bool output_csv)
{
	// Output markdown file
	std::fstream output_file {path(ROOT_DIR) / path("outputs/benchmark.md"), std::ios::out | std::ios::trunc};
	if(!output_file) {
		std::cerr << "Canont open \"outputs/benchmark.md\"" << std::endl;
		return;
	}

	// Creating one thread per benchmark
	std::vector<std::thread> workers;

	workers.emplace_back( [&](){
		auto gaussian_square = load_test_2d("gaussian_square");
		auto gaussian_square_res = benchmark_all_opt(gaussian_square, gen_lambda_prop(-5.0, -1.0, 0.5));
		print_benchmark_md(gaussian_square_res, output_file);
		if(output_csv) save_log_densities(gaussian_square_res);
	});

	workers.emplace_back([&](){	
		auto infections_southampton = load_test_2d("infections_southampton");
		auto infections_southampton_res = benchmark_all_opt(infections_southampton, gen_lambda_prop(-5.0, -1.0, 1.0));
		print_benchmark_md(infections_southampton_res, output_file);
		if(output_csv) save_log_densities(infections_southampton_res);
	});

	workers.emplace_back([&](){
		auto horseshoe = load_test_2d("horseshoe");
		auto horseshoe_res = benchmark_all_opt(horseshoe, gen_lambda_prop(-5.0, -1.0, 1.0));
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

	workers.emplace_back([&](){
		auto snp500 = load_test_snp500("snp500");
		auto snp500_res = benchmark_all_opt(snp500, gen_lambda_prop(-2.0, 2.0, 1.0));
		print_benchmark_md(snp500_res, output_file);
		if(output_csv) save_log_densities(snp500_res);
	});

	for(auto &t: workers) {
		t.join();
	}
}

} // namespace de