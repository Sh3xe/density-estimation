#include "space_time_params.hpp"
#include "utilities.hpp"
#include <fdaPDE/fdapde.h>

using namespace fdapde;
namespace fs = std::filesystem;
using namespace Eigen;

constexpr double LINK_TOL = 1e-3;
constexpr int MAX_ITER = 100;
constexpr int CV_K = 5;
constexpr double STEP_SIZE = 1e-2;
constexpr double ERR_TOL = LINK_TOL;

struct TestScenario {
	Triangulation<2, 2> space_domain;
	Triangulation<1, 1> time_domain;
	MatrixXd g_init;
	MatrixXd dataset;
};

struct EvaluationCounter {
	int func_eval_count = 0;
	int grad_eval_count = 0;

	void reset() {
		func_eval_count = 0;
		grad_eval_count = 0;
	}

	template <typename Opt, typename Obj> bool grad_hook(Opt& opt, Obj& obj) {
		++grad_eval_count;
		return false;
	}

	template <typename Opt, typename Obj> bool eval_hook(Opt& opt, Obj& obj) {
		++func_eval_count;
		return false;
	}
};

TestScenario load_test() {
	// Load the test scenario
	auto dir = path(ROOT_DIR) / path("data/unit_square_space_time");
	Triangulation<2, 2> space_domain(
		(dir / path("mesh_vertices.csv")).string(),
		(dir / path("mesh_elements.csv")).string(),
		(dir / path("mesh_boundary.csv")).string(),
		true, true
	);

	Triangulation<1, 1> time_domain = Triangulation<1, 1>::UnitInterval(7);

	MatrixXd g_init = read_csv<double>((dir / path("f_init.csv")).string())
		.as_matrix()
		.array()
		.log();

	MatrixXd dataset(500, 3);
	dataset.leftCols(2)  = read_csv<double>((dir / path("sample_space.csv")).string()).as_matrix();
	dataset.rightCols(1) = read_csv<double>((dir / path("sample_time.csv")).string()).as_matrix();

	TestScenario scenario;
	scenario.space_domain = space_domain;
	scenario.time_domain = time_domain;
	scenario.g_init = g_init;
	scenario.dataset = dataset;

	std::cout << space_domain.nodes().rows() << " " << space_domain.nodes().cols() << std::endl;
	std::cout << time_domain.nodes().rows() << " " << time_domain.nodes().cols() << std::endl;

	return scenario;
}

template < typename DEPDESolver >
double cv_error(
	const double lambda_space,
	const double lambda_time,
	DEPDESolver &solver, // fe_de_separable
	TestScenario &scenario,
	int cv_k
) {
	double err_tot = 0.0;
	
	for(int i = 0; i < cv_k; ++i) {
		// Split the dataset & setup solver
		const auto [train_set, test_set] = utils::split_dataset(scenario.dataset, CV_K, i);
		GeoFrame geo_data(scenario.space_domain, scenario.time_domain);

		auto& data_layer = geo_data.template insert_scalar_layer<POINT, POINT>("l1", train_set);
		DEPDE<typename DEPDESolver::solver_t> model(geo_data, solver);

		model.set_llik_tolerance(LINK_TOL);

		// Solve
		model.fit(
			lambda_space, lambda_time, scenario.g_init,
			LBFGS<fdapde::Dynamic>(MAX_ITER, ERR_TOL, STEP_SIZE, 10), 
			WolfeLineSearch()
		);

		// -------- SPACE-ONLY ERROR COMPUTATION --------
		// Compute the error on the test set
		// FeSpace Vh(scenario.discretization, P1<1>);
		// FeFunction log_dens_func(Vh, model.log_density());

		// double total_dens_eval = 0.0;
		// for(int i = 0; i < test_set.rows(); ++i) {
		// 	total_dens_eval += std::exp( log_dens_func(test_set.row(i)) );
		// }

		// err_tot += integral(scenario.discretization, QS2DP4)(exp(2.0*log_dens_func));
		// err_tot	-= 2.0 * (total_dens_eval / (double)test_set.rows());
	}

	return err_tot / static_cast<double>(CV_K);
}

double fit_space_time(TestScenario &scenario, const double lambda_space, const double lambda_time) {
	// Defining the model
	FeSpace Vh(scenario.space_domain, P1<1>);   // linear finite element in space
	TrialFunction f(Vh);
	TestFunction  v(Vh);
	auto a_D = integral(scenario.space_domain)(dot(grad(f), grad(v)));
	ZeroField<2> u_D;
	auto F_D = integral(scenario.space_domain)(u_D * v);

	BsSpace Qh(scenario.time_domain, 3);   // cubic B-spline in time
	TrialFunction g(Qh);
	TestFunction  w(Qh);
	auto a_T = integral(scenario.time_domain)(dxx(g) * dxx(w));
	ZeroField<1> u_T;
	auto F_T = integral(scenario.time_domain)(u_T * w);

	// Model
	auto solver = fe_de_separable(std::pair {a_D, F_D}, std::pair {a_T, F_T});
	return cv_error(lambda_space, lambda_time, solver, scenario, CV_K);
}

void test() {
	TestScenario scenario = load_test();
	fit_space_time(scenario, 0.00025, 0.01);
}