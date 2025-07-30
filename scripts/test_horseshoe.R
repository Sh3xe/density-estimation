library(fdaPDE)
library(viridis)
rm(list = ls())

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")

# Define the domain
nodes <- boundary.domain.2()
# Remove duplicate coordinates
nodes <- nodes[-nrow(nodes),]
nodes <- nodes[-45,]

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
plot.mesh.2D(mesh = mesh)

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
	lambda <- lambda_proposal,
	nfolds = 10,
	tol1 = 1e-4,
	nsimulations = 500,
	preprocess_method = "RightCV",
	step_method = "Wolfe_Method",
	direction_method = "L-BFGS10"
)

best_lambda_id <- match(de$lambda, lambda_proposal)

write.csv(de$f_init[,best_lambda_id], file.path(output_dir, "f_init.csv"))

# Calculate the true density function
true_density <- dens.func.1(mesh$nodes)

write.csv(true_density, file.path(output_dir, "true_density.csv"))

# Plot
X <- seq(from = min(boundary.domain.2()[,1]), to = max(boundary.domain.2()[,1]),
         length.out = 100)
Y <- seq(from = min(boundary.domain.2()[,2]), to = max(boundary.domain.2()[,2]),
         length.out = 100)

grid <- expand.grid(X, Y)
mesh.eval <- create.mesh.2D(nodes = grid)

# Set up the finite element basis
FEMbasis.eval <- create.FEM.basis(mesh = mesh.eval)

estimated_log_dens <- read.csv("outputs/cpp_lbfgs30_horseshoe_log_density.csv")
FEMfunction <- FEM(coeff = estimated_log_dens, FEMbasis = FEMbasis)
evaluation <- eval.FEM(FEM = FEMfunction, locations = mesh.eval$nodes)

# Compute the density estimate on the nodes of the finer mesh
estimated_density <- exp(evaluation)

# True density
true_density <- dens.func.2(data = mesh.eval$nodes,
                            mesh = mesh.eval)  # function in helper_functions_data.R

# Plot
par(mfrow = c(1,2), mai = c(0.5,0.25,0.5,0.5))

# function in helper_functions_plot.R
plot.density.2D(X = X, Y = Y, Z = true_density, boundary = boundary.domain.2(), colorscale="viridis")

# function in helper_functions_plot.R
plot.density.2D(X = X, Y = Y, Z = estimated_density, boundary = boundary.domain.2(), colorscale="viridis")
