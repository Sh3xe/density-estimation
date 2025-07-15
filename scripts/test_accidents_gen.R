library(fdaPDE)
rm(list = ls())
directory <- "data/accidents_bergamo"

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")

# Define the domain
load(file = file.path(directory, "domain.RData"))
mesh$nodesmarkers[mesh$nodesmarkers == FALSE] <- 0
mesh$nodesmarkers[mesh$nodesmarkers == TRUE] <- 1
submesh$nodesmarkers[submesh$nodesmarkers == FALSE] <- 0
submesh$nodesmarkers[submesh$nodesmarkers == TRUE] <- 1

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = submesh)

# Generate the data
load(file = file.path(directory,"data.RData"))
data = data[,1:2]

# Save the mesh and samples
write.csv(mesh$nodes, file.path(directory,"mesh_vertices.csv"))
write.csv(mesh$nodesmarker, file.path(directory, "mesh_boundary.csv"))
write.csv(mesh$edges, file.path(directory, "mesh_edges.csv"))
write.csv(data, file.path(directory, "sample.csv"))

# Calculate f_init
lambda_proposal <- 10^seq(from = -1, to = -5, by = -0.5)
de <- DE.FEM(
	data = data,
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
