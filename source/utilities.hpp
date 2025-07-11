#pragma once

#include <Eigen/Dense>
#include <chrono>
#include <fdaPDE/fdapde.h>

using matrix_t = Eigen::Matrix<double, fdapde::Dynamic, fdapde::Dynamic>;
using vector_t = Eigen::Matrix<double, fdapde::Dynamic, 1>;
using microsec = std::chrono::duration<double, std::micro>;

std::string duration_to_str(const microsec &duration);