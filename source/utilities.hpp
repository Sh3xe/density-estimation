#pragma once

#include <chrono>

using microsec = std::chrono::duration<double, std::micro>;

std::string duration_to_str(const microsec &duration);