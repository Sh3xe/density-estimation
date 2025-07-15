library(fdaPDE)
rm(list = ls())

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")

# Define the domain
nodes <- boundary.domain.2()
segments <- matrix(
	data = c(1, rep(2:nrow(nodes), each = 2), 1),
  nrow = nrow(nodes),
	ncol = 2,
	byrow = TRUE
)

# Create a regular mesh of the horseshoe-shaped domain
mesh <- create.mesh.2D(nodes = nodes, segments = segments, order = 1)
mesh <- refine.mesh.2D(mesh = mesh, maximum_area = 0.025)

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Generate the data
n <- 1000
data <- generate.spatial.data.2(N = n, mesh = mesh)     # function in helper_functions_data.R

# Save the mesh and samples
output_dir <- "data/horseshoe"

write.csv(mesh$nodes, file.path(output_dir,"mesh_vertices.csv"))
write.csv(mesh$triangles, file.path(output_dir, "mesh_elements.csv"))
write.csv(mesh$nodesmarker, file.path(output_dir, "mesh_boundary.csv"))
write.csv(data, file.path(output_dir, "sample.csv"))

# Calculate f_init
lambda_proposal <- 10^seq(from = -1, to = -5, by = -1)

de <- DE.FEM(
	data = data,
	FEMbasis = FEMbasis,
	lambda <- 0.01,
	nfolds = 10,
	tol1 = 1e-4,
	nsimulations = 500,
	# preprocess_method = "RightCV",
	step_method = "Wolfe_Method",
	direction_method = "L-BFGS10"
)

best_lambda_id <- match(de$lambda, lambda_proposal)

write.csv(de$f_init, file.path(output_dir, "f_init.csv"))

# Calculate the true density function
true_density <- dens.func.1(mesh$nodes)

write.csv(true_density, file.path(output_dir, "true_density.csv"))