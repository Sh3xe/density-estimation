#pragma once

#include <Eigen/Dense>
#include <chrono>

using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using microsec = std::chrono::duration<double, std::micro>;

std::string duration_to_str(const microsec &duration);