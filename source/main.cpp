#include "benchmark_de.hpp"
#include "benchmark_low_dim.hpp"
#include <cstring>

int main(int argc, char **argv) {
	bool output_csv = (argc > 1 && strcmp(argv[1], "output_csv") == 0);
	std::string type = "de";
	if(argc > 2) {
		type = argv[2];
	}

	if(type == "de") {
		std::cout << "Starting DE benchmarks, \"output_csv\" = " << (output_csv ? "Yes": "No") << std::endl;
		de::full_benchmark(output_csv);
	} else if( type == "lowdim") {
		std::cout << "Starting low dim benchmarks, \"output_csv\" = " << (output_csv ? "Yes": "No") << std::endl;
		lowdim::full_benchmark(output_csv);
	} else {
		std::cout << "No type specified, type must be \"lowdim\" or \"de\"" << std::endl;
	}

	return 0;
}