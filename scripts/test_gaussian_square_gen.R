library(fdaPDE)
rm(list = ls())

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")

# Define the domain
X <- seq(from = 0, to = 1, length.out = 11)
Y <- seq(from = 0, to = 1, length.out = 11)
grid <- expand.grid(X, Y)
bounds <- grid[(grid$Var1 %in% c(0, 1)) | (grid$Var2 %in% c(0, 1)), ]

# Create a regular mesh of the unit square
mesh <- create.mesh.2D(nodes = bounds, order = 1)
mesh <- refine.mesh.2D(mesh = mesh, maximum_area = 0.0025, minimum_angle = 25)

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Generate the test data
N <- 1000
gaussian_data <- generate.spatial.data.1(N)

# Save the mesh and samples
directory <- "data/gaussian_square"

write.csv(mesh$nodes, file.path(directory,"mesh_vertices.csv"))
write.csv(mesh$triangles, file.path(directory, "mesh_elements.csv"))
write.csv(mesh$nodesmarker, file.path(directory, "mesh_boundary.csv"))
write.csv(gaussian_data, file.path(directory, "sample.csv"))

# Calculate f_init
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

best_lambda_id <- match(de$lambda, lambda_proposal)

write.csv(de$f_init[,best_lambda_id], file.path(directory, "f_init.csv"))

# Calculate the true density function
true_density <- dens.func.1(mesh$nodes)

write.csv(true_density, file.path(directory, "true_density.csv"))