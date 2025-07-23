#pragma once

#include "cmake_defines.hpp"

#include <string>
#include <chrono>
#include <filesystem>
#include <vector>
#include <Eigen/Dense>

using path = std::filesystem::path;
using microsec = std::chrono::duration<double, std::micro>;

namespace utils {
	
	Eigen::MatrixXd load_csv(const std::string &filepath);
	void write_csv(const std::string &filepath, const std::vector<std::string> &columns, const Eigen::MatrixXd &matrix);
	std::string duration_to_str(const microsec &duration);
	double mean(const Eigen::MatrixXd &df, size_t col_idx);
	double std(const Eigen::MatrixXd &df, size_t col_idx);

} // namespace utils