library(fdaPDE)
library(viridis)
rm(list = ls())
directory <- "data/infections_southampton"

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")

# Load data
load(file = file.path(directory, "infections_southampton.RData"))

# Define the domain
nodes <- cbind(as.data.frame(xyt$window$bdry[[1]]$x), as.data.frame(xyt$window$bdry[[1]]$y))
segments <- cbind(1:dim(nodes)[1], c(2:dim(nodes)[1], 1))

# Define the mesh
mesh <- create.mesh.2D(nodes = nodes, segments = segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 4.25, minimum_angle = 25)

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Import the data
data <- cbind(xyt$x, xyt$y)

# Save the mesh and samples
write.csv(data, file.path(directory, "sample.csv"))
write.csv(mesh$nodes, file.path(directory,"mesh_vertices.csv"))
write.csv(mesh$triangles, file.path(directory, "mesh_elements.csv"))
write.csv(mesh$nodesmarker, file.path(directory, "mesh_boundary.csv"))

# Calculate f_init
lambda_proposal <- 10^seq(from = -1, to = -5, by = -0.5)

de <- DE.FEM(
	data=data,
	FEMbasis = FEMbasis,
	lambda = lambda_proposal,
	nfolds = 10,
	tol1 = 1e-4,
	nsimulations = 1000,
	preprocess_method = "RightCV",
	step_method = "Wolfe_Method",
	direction_method = "L-BFGS10"
)

best_lambda_id <- match(de$lambda, lambda_proposal)

write.csv(de$f_init[,best_lambda_id], file.path(directory, "f_init.csv"))

# Load the data
f_init <- read.csv(file.path(directory, "f_init.csv"))
sample <- read.csv(file.path(directory, "sample.csv"))
log_dens <- read.csv("./outputs/cpp_lbfgs30_infections_southampton_log_density.csv")

# Plot
FEMfunction <- FEM(coeff = log_dens$V0, FEMbasis = FEMbasis)
evaluation <- eval.FEM(FEM = FEMfunction, locations = mesh$nodes)
estimated_density <- exp(evaluation)
plot.density.2D.map(coeff = estimated_density, mesh = mesh, colorscale = "viridis")

