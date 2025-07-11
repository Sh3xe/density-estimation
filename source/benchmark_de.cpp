#include "benchmark_de.hpp"
#include <filesystem>

using path = std::filesystem::path;
using namespace fdapde;

std::vector<DEBenchmarkResult> benchmark_de_problem(
	const std::string &mesh_dir_str,
	const std::string &title
) {
	path dir = path(mesh_dir_str);

	// Geometry
	Triangulation<2, 2> domain(
		(dir/path("points.csv")).string(),
		(dir/path("elements.csv")).string(),
		(dir/path("boundary.csv")).string(),
		true, true
	);

	// Data
	vector_t g_init = read_csv<double>((dir/path("f_init.csv")).string())
		.as_matrix()
		.array()
		.log();

	GeoFrame data(domain);
	auto& l1 = data.insert_scalar_layer<POINT>("l1", (dir/path("points.csv")).string());

	// Physics
	FeSpace Vh(domain, P1<1>);
	TrialFunction f(Vh);
	TestFunction  v(Vh);
	auto a = integral(domain)(dot(grad(f), grad(v)));
	ZeroField<2> u;
	auto F = integral(domain)(u * v);

	// Modeling
	DEPDE solver(data, fe_de_elliptic(a, F));
	solver.set_llik_tolerance(1e-8);

	// Benchmark
	std::vector<DEBenchmarkResult> res;
	res.reserve(9);
	
	auto lbfgs30 = LBFGS<Dynamic, WolfeLineSearch> {500, 1e-8, 1e-2, 30};
	res.push_back(benchmark_solver(solver, lbfgs30, g_init, title));

	auto cg_pr = ConjugateGradient<Dynamic, WolfeLineSearch> {500, 1e-5, 1e-2, true};
	res.push_back(benchmark_solver(solver, cg_pr, g_init, title));

	return res;
}