#include "space_time_params.hpp"
#include "utilities.hpp"
#include <vector>
#include <fdaPDE/optimization.h>

using namespace fdapde;

double bilinear_interpolation(
	const Eigen::MatrixXd& grid,
	const Eigen::VectorXd& x_grid,
	const Eigen::VectorXd& y_grid,
	double x, double y) {

	// Find the indices of the four neighboring points
	int x1_idx = 0, x2_idx = 0;
	int y1_idx = 0, y2_idx = 0;

	// Find x indices
	for (int i = 0; i < x_grid.size() - 1; ++i) {
		if (x_grid(i) <= x && x_grid(i+1) > x) {
			x1_idx = i;
			x2_idx = i + 1;
		}
	}

	// Find y indices
	for (int i = 0; i < y_grid.size() - 1; ++i) {
		if (y_grid(i) <= y && y_grid(i+1) > y) {
			y1_idx = i;
			y2_idx = i + 1;
		}
	}

	// Get the four nearest values
	double Q11 = grid(x1_idx, y1_idx);
	double Q12 = grid(x1_idx, y2_idx);
	double Q21 = grid(x2_idx, y1_idx);
	double Q22 = grid(x2_idx, y2_idx);

	// Calculate interpolation weights
	double x_frac = (x - x_grid(x1_idx)) / (x_grid(x2_idx) - x_grid(x1_idx));
	double y_frac = (y - y_grid(y1_idx)) / (y_grid(y2_idx) - y_grid(y1_idx));

	// Perform bilinear interpolation
	double interpolated_value = 
		(1 - x_frac) * (1 - y_frac) * Q11 +
		(1 - x_frac) * y_frac * Q12 +
		x_frac * (1 - y_frac) * Q21 +
		x_frac * y_frac * Q22;

	return interpolated_value;
}

struct ParamFunction {
	int max_evals = 100;
	bool log = false;
	Eigen::VectorXd space_values;
	Eigen::VectorXd time_values;
	Eigen::VectorXd cv_values;
	Eigen::MatrixXd cv_matrix;
	double x_min = 0.0, x_max = 1.0;
	double y_min = 0.0, y_max = 1.0;
	
	int eval_counter = 0;
	double best_value = std::numeric_limits<double>::max();
	std::vector<double> values;
	std::vector<Eigen::VectorXd> eval_points;

	ParamFunction(bool log = true): log(log) {
		auto dir = path(ROOT_DIR);
		space_values = utils::load_csv(dir / path("outputs/st_hm_space.csv"), 0, 0);
		time_values = utils::load_csv(dir / path("outputs/st_hm_time.csv"), 0, 0);
		cv_values = utils::load_csv(dir / path("outputs/st_hm_cverr.csv"), 0, 0);

		cv_matrix.resize(space_values.size(), time_values.size());
		for (int i = 0; i < space_values.size(); ++i) {
			for (int j = 0; j < time_values.size(); ++j) {
				cv_matrix(i, j) = cv_values[i * time_values.size() + j];
			}
		}

		x_min = space_values.minCoeff(); x_max = space_values.maxCoeff();
		y_min = time_values.minCoeff(); y_max = time_values.maxCoeff();
	}

	double operator()(const Eigen::VectorXd &input) {

		double x = std::pow(10, input(0));
		double y = std::pow(10, input(1));
		
		double x_diff = 0.0;
		double y_diff = 0.0;

		if(x < x_min || x > x_max) {
			x_diff = std::min(std::abs(x-x_min), std::abs(x-x_max));
			x = std::max(std::min(x_max - 0.1, x), x_min + 0.1);
		}
		if(y < y_min || y > y_max) {
			y_diff = std::min(std::abs(y-y_min), std::abs(y-y_max));
			y = std::max(std::min(y_max - 0.1, y), y_min + 0.1);
		}

		double current_eval = bilinear_interpolation(cv_matrix, space_values, time_values, x, y);
		current_eval += x_diff * x_diff * 10;
		current_eval += y_diff * y_diff * 10;

		if(log) {
			std::cout << "Func(" << input(0) << "," << input(1) << ") = " << current_eval << std::endl;
		}

		eval_counter++;
		values.push_back(current_eval);
		eval_points.push_back(input);

		if(best_value > current_eval) {
			best_value = current_eval;
		}

		return current_eval;
	}

	void save_data(const std::string &file_name) {
		Eigen::MatrixXd df = Eigen::MatrixXd::Zero(values.size(), 3);
		for(int i = 0; i < values.size(); ++i) {
			df(i, 0) = values[i];
			df(i, 1) = eval_points[i](0);
			df(i, 2) = eval_points[i](1);
		}
		utils::write_csv(file_name, {"value", "x", "y"}, df);
	}

	template <typename Opt>
	bool stop_if(Opt& optimizer) {
		if(log) {
			std::cout << "stop_if call " << std::endl;
		}
		return eval_counter >= max_evals;
	}

	void reset() {
		eval_counter = 0;
		values.clear();
		eval_points.clear();
		best_value = std::numeric_limits<double>::max();
	}
};

void test() {
	ParamFunction func {false};
	Eigen::VectorXd x0(2);
	x0(0) = log10(func.space_values(func.space_values.size() / 2));
	x0(1) = log10(func.time_values(func.time_values.size() / 2));
	std::cout << "x0 = " << x0(0) << ", " << x0(1) << std::endl;

	{
		GeneticOptim<fdapde::Dynamic> go{100, 0.0, 100, 5};
		func.reset();
		go.optimize(func, x0, GaussianMutation(0.1, 0.95), BinaryTournamentSelection());
		std::cout << "Min: " <<  go.value() << ", Func evals: " << func.eval_counter << std::endl;
		func.save_data("outputs/cal_go_gm_bts.csv");
	}

	{
		GeneticOptim<fdapde::Dynamic> go{100, 0.0, 100, 5};
		func.reset();
		go.optimize(func, x0, GaussianMutation(0.1, 0.95), CrossoverMutation(), BinaryTournamentSelection());
		std::cout << "Min: " <<  go.value() << ", Func evals: " << func.eval_counter << std::endl;
		func.save_data("outputs/cal_go_gm_cm_bts.csv");
	}

	{
		GeneticOptim<fdapde::Dynamic> go{100, 0.0, 100, 5};
		func.reset();
		go.optimize(func, x0, GaussianMutation(0.1, 0.95), RankSelection());
		std::cout << "Min: " <<  go.value() << ", Func evals: " << func.eval_counter << std::endl;
		func.save_data("outputs/cal_go_gm_rs.csv");
	}

	{
		GeneticOptim<fdapde::Dynamic> go{100, 0.0, 100, 5};
		func.reset();
		go.optimize(func, x0, GaussianMutation(0.1, 0.95), CrossoverMutation(), RankSelection());
		std::cout << "Min: " <<  go.value() << ", Func evals: " << func.eval_counter << std::endl;
		func.save_data("outputs/cal_go_gm_cm_rs.csv");
	}

	{
		NelderMead<fdapde::Dynamic> nm {100, 0.0};
		func.reset();
		nm.optimize(func, x0, GaussianMutation(0.1, 0.95), BinaryTournamentSelection());
		std::cout << "Min: " <<  nm.value() << ", Func evals: " << func.eval_counter << std::endl;
		func.save_data("outputs/cal_nm.csv");
	}
}