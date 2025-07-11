#include "utilities.hpp"
#include <string>

std::string duration_to_str(const microsec &duration) {
	if(duration.count() < 1e3) {
		return std::to_string((int)duration.count()) + std::string("\\textmu s");
	} else if(duration.count() < 1e6) {
		return std::to_string((int)(duration.count()*1e-3) ) + std::string("ms");
	} else if(duration.count() < 1e9) {
		return std::to_string((int)(duration.count()*1e-6) ) + std::string("s");
	} else {
		return ">1000s";
	}
}