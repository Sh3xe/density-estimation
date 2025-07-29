#include "benchmark_de.hpp"
#include <cstring>

int main(int argc, char **argv) {
	bool output_csv = (argc > 1 && strcmp(argv[1], "output_csv") == 0);
	std::cout << "Starting DE benchmarks, \"output_csv\" = " << (output_csv ? "Yes": "No") << std::endl;
	de::full_benchmark(output_csv);
	return 0;
}