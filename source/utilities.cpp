#include "utilities.hpp"
#include <cmath>
#include "rapidcsv.h"
#include <fstream>

std::string utils::duration_to_str(const microsec &duration) {
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

Eigen::MatrixXd utils::load_csv(const std::string &filepath) {
	rapidcsv::Document doc((path(ROOT_DIR) / path(filepath)).string(), rapidcsv::LabelParams(-1, -1));
	Eigen::MatrixXd csv_mat = Eigen::MatrixXd::Zero(doc.GetRowCount(), doc.GetColumnCount());

	for(int i = 0; i < csv_mat.rows(); ++i) {
		auto vec = doc.GetRow<double>(i);
		for(int j = 0; j < csv_mat.cols(); ++j)
			csv_mat(i, j) = vec[j];
	}

	return csv_mat;
}

void utils::write_csv(
	const std::string &filepath,
	const std::vector<std::string> &columns,
	const Eigen::MatrixXd &matrix) {
	assert(columns.size() == matrix.cols());
	std::ofstream file(  path(ROOT_DIR) / path(filepath), std::ios::out | std::ios::trunc );
	if(!file) {
		std::cerr << "Cannot save csv file to " << filepath << std::endl;
		return;
	}
	
	// Header
	for(const auto &col_name: columns)
	file << "," << col_name;
	file << "\n";
	
	// content
	for(size_t row_idx = 0; row_idx < matrix.rows(); ++row_idx) {
		file << row_idx;
		for(size_t col_idx = 0; col_idx < columns.size(); ++col_idx) {
			file << "," << matrix(row_idx, col_idx);
		}
		file << "\n";
	}

	file << std::flush;
}

double utils::mean(const Eigen::MatrixXd &df, size_t col_idx) {
	double tot = 0.0;
	for(size_t j = 0; j < df.rows(); ++j) {
		tot += df(j, col_idx);
	}
	return tot / (double)df.rows();
}

double utils::std(const Eigen::MatrixXd &df, size_t col_idx) {
	double df_mean = utils::mean(df, col_idx);
	double tot = 0.0;

	for(size_t j = 0; j < df.rows(); ++j) {
		double x = df(j,col_idx) - df_mean;
		tot += x*x;
	}

	return std::sqrt( tot / (double)(df.rows()-1) );
}