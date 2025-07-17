library(fdaPDE)
rm(list = ls())
directory <- "data/curved"

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")

# Define the domain
vertices <- read.table(file.path(directory, "mesh_vertices.txt"), quote = "\"", comment.char = "")
triangles <- read.table(file.path(directory, "mesh_triangles.txt"), quote = "\"", comment.char = "")[,1:3]

# Create the mesh
mesh <- create.mesh.2.5D(nodes = vertices[,1:3], triangles = triangles[,1:3])

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Generate the data
n <- 1000
data <- generate.spatial.data.4(N = n, mesh = mesh)

# Save the mesh and samples
write.csv(mesh$nodes, file.path(directory,"mesh_vertices.csv"))
write.csv(mesh$triangles, file.path(directory, "mesh_elements.csv"))
write.csv(mesh$nodesmarker, file.path(directory, "mesh_boundary.csv"))
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

# Calculate the true density function
true_density <- dens.func.3(data = mesh$nodes)
write.csv(true_density, file.path(directory, "true_density.csv"))