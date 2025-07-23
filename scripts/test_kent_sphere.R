library(fdaPDE)
rm(list = ls())
directory <- "data/kent_sphere"

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")

# Define the domain
vertices <- read.table(file.path(directory, "mesh_vertices.txt" ), quote = "\"", comment.char = "")
triangles <- read.table(file.path(directory, "mesh_triangles.txt" ), quote = "\"", comment.char = "")

# Create a regular mesh of the unit square
mesh <- create.mesh.2.5D(nodes = vertices[,1:3], triangles = triangles[,1:3])
# mesh$nodesmarker[mesh$nodesmarker == TRUE] = 1
# mesh$nodesmarker[mesh$nodesmarker == FALSE] = 0

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Generate the test data
n <- 1000
data <- generate.spatial.data.3(N = n)

# Save the mesh and samples
directory <- "data/kent_sphere"

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

# Ploting
estimated_density <- read.csv("outputs/lbfgs30_kent_sphere_log_density.csv")
FEMfunction <- FEM(coeff = estimated_density$V0, FEMbasis = FEMbasis)
evaluation <- eval.FEM(FEM = FEMfunction, locations = mesh$nodes)

# Compute the density estimate on the nodes of the finer mesh
estimated_density <- exp(evaluation)

# True density
true_density <- dens.func.3(data = mesh.eval$nodes) # function in helper_functions_data.R

m = min(true_density, estimated_density, na.rm = TRUE)
M = max(true_density, estimated_density, na.rm = TRUE)

# Plot
plot.density.2.5D(coeff = true_density, mesh = mesh, min_range = m, max_range = M,
                  colorscale = "viridis")     # function in helper_functions_plot.R
