#include "benchmark_de.hpp"
#include "benchmark_low_dim.hpp"

#include <cstring>
#include <string>
#include <iostream>

int main(int argc, char **argv) {
	std::string type = "";
	if(argc > 1) {
		type = argv[1];
	}

	bool output_csv = (argc > 2 && strcmp(argv[2], "output_csv") == 0);

	if(type == "de") {
		std::cout << "Starting DE benchmarks, \"output_csv\" = " << (output_csv ? "Yes": "No") << std::endl;
		de_full_benchmark(output_csv);
	} else if( type == "lowdim") {
		std::cout << "Starting low dim benchmarks, \"output_csv\" = " << (output_csv ? "Yes": "No") << std::endl;
		lowdim_full_benchmark(output_csv);
	} else if( type == "compde") {
		print_l2_errors("gaussian_square");
	} else {
		std::cout << "No type specified, type must be \"lowdim\" or \"de\": [name] [type] output_csv?" << std::endl;
	}

	return 0;
}