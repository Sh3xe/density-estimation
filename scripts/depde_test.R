# Import the library
library(fdaPDE)           # v. 1.1-21 (2025)
rm(list = ls())

setwd("~/Documents/Boulot/PRE/density-estimation/")

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")


## [DOMAIN]
# Define the bounds of the unit square
X <- seq(from = 0, to = 1, length.out = 11)
Y <- seq(from = 0, to = 1, length.out = 11)
grid <- expand.grid(X, Y)
bounds <- grid[(grid$Var1 %in% c(0, 1)) | (grid$Var2 %in% c(0, 1)), ]

## [MESH]
# Create a regular mesh of the unit square
mesh <- create.mesh.2D(nodes = bounds, order = 1)
mesh <- refine.mesh.2D(mesh = mesh, maximum_area = 0.0025, minimum_angle = 25)

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Plot
plot.mesh.2D(mesh = mesh)   # function in helper_functions_plot.R

## [SIMPLE GAUSSIAN DATA]
N <- 1000
gaussian_data <- generate.spatial.data.1(N)

plot.sample.2D(mesh=mesh, data=gaussian_data)

# calculate 
output_dit <- "data/gaussian_2d"
lambda_proposal <- 10^seq(from = -1, to = -5, by = -0.5)


de <- DE.FEM(
	data = gaussian_data,
	FEMbasis = FEMbasis,
	lambda <- lambda_proposal,
	nfolds = 10,
	tol1 = 1e-4,
	nsimulations = 1000,
	preprocess_method = "RightCV",
	step_method = "Wolfe_Method",
	direction_method = "L-BFGS10"
)

