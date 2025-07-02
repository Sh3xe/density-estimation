cmake_minimum_required(VERSION 3.14)

project(advanced_optim_applic)

add_executable( advanced_optim_applic
	source/main.cpp
)

target_compile_features( advanced_optim_applic PUBLIC cxx_std_20)

target_include_directories(advanced_optim_applic PUBLIC 
	fdaPDE-cpp
	fdaPDE-cpp/fdaPDE/core
)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
 
target_link_libraries(advanced_optim_applic Eigen3::Eigen)