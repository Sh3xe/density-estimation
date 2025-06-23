rm(list=ls())
source("read.mesh.R")
library(fdaPDE)

# Loading Data
mesh_path <- "./data/square_16.mesh"
point_path <- "./data/test_simple/points.csv"
domain <- read.mesh(mesh_path)
domain$nodes <- domain$nodes * 6 - 3
mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements, 
                       segments = domain$faces)
data <- read.csv(point_path)
data <- data[, c("data_x", "data_y")]

# Model fiting
de_fem_basis <- create.FEM.basis(mesh)

# Importing
result_path <- "./data/test_simple/g_cpp_lbfgs30_l0.1_wolfe.csv"
g <- read.csv(result_path)
g <- g$x

# Evaluation
evals_pt <- eval.FEM(FEM(de_fem_basis, coeff = g), locations = data)
error <- -sum(evals_pt)
sprintf("%.10f",error)
