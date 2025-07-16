#include "benchmark_de.hpp"

void print_benchmark_md( const std::vector<DEBenchmarkResult> &benchmark, std::ostream &os ) {
	os << "### " << benchmark[0].test_title << "\n";
	os << "Optimizer|Iterations|Loss|CV error|lambda|\n";
	os << "-|-|-|-|-\n";

	for(const auto &res: benchmark) {
		os << res.optimizer << "|" << res.iterations << "|" << res.loss << "|" << res.cv_error << "|" << res.lambda << "\n";
	}
}

void full_benchmark()
{
	// Load test scenario
	// auto snp500 = load_test_snp500("snp500");
	auto gaussian_square = load_test_scenario("gaussian_square");
	auto infections_southampton = load_test_scenario("infections_southampton");
	// auto horseshoe = load_test_scenario("horseshoe");
	// auto kent_sphere = load_test_scenario("kent_sphere");
	// auto gaussian_curved = load_test_scenario("gaussian_curved");
	// auto earthquake = load_test_scenario("earthquake");
	// auto accidents_bergamo = load_test_scenario("accidents_bergamo");

	// Optimizers
	auto lbfgs30 = std::make_pair(std::string{"lbfgs30"}, LBFGS<Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, 30, ERR_TOL, STEP_SIZE});
	auto grad_descent = std::make_pair(std::string{"gd"}, GradientDescent<Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE});
	// auto cg_pr = std::make_pair( std::string{""}, ConjugateGradient<Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, true});
	// auto cg_fr = std::make_pair( std::string{""}, ConjugateGradient<Dynamic, WolfeLineSearchNew> {MAX_ITERATIONS, ERR_TOL, STEP_SIZE, false});
	// auto nelder_mead = std::make_pair( std::string{""}, NelderMead<Dynamic> {ERR_TOL, STEP_SIZE});

	// DE problem formulation
	FeSpace Vh(space_discr, P1<1>);
	TrialFunction f(Vh);
	TestFunction  v(Vh);
	auto a = integral(space_discr)(dot(grad(f), grad(v)));
	ZeroField<2> u;
	auto F = integral(space_discr)(u * v);
	auto solver = fe_de_elliptic(a, F);

	// Benchmarks
	auto gaussian_square_res = benchmark_all(gaussian_square, solver, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	print_benchmark_md(gaussian_square_res);

	auto infections_southampton_res = benchmark_all(infections_southampton, solver, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	print_benchmark_md(infections_southampton_res);

	// auto horseshoe_res = benchmark_all(horseshoe, solver, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(horseshoe_res);

	// auto kent_sphere_res = benchmark_all(kent_sphere, solver, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(kent_sphere_res);

	// auto gaussian_curved_res = benchmark_all(gaussian_curved, solver, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(gaussian_curved_res);

	// auto earthquake_res = benchmark_all(earthquake, solver, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(earthquake_res);

	// auto accidents_bergamo_res = benchmark_all(accidents_bergamo, solver, lbfgs30, grad_descent, /*cg_pr, cg_fr, nelder_mead*/);
	// print_benchmark_md(accidents_bergamo_res);

	return 0;
}