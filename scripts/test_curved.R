# Import the library
library(fdaPDE)
rm(list = ls())
directory <- "data/curved"

source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")

# Import vertices and triangles
vertices <- read.table(file.path(directory, "mesh_vertices.txt"), quote = "\"", comment.char = "")
triangles <- read.table(file.path(directory, "mesh_triangles.txt"), quote = "\"", comment.char = "")[,1:3]
mesh$nodesmarker[mesh$nodesmarker == TRUE] = 1
mesh$nodesmarker[mesh$nodesmarker == FALSE] = 0

# Create a regular mesh of the surface of the unit sphere
mesh <- create.mesh.2.5D(nodes = vertices[,1:3], triangles = triangles[,1:3])

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Sample size
n <- 1000

# Generate the data
data <- generate.spatial.data.4(N = n, mesh = mesh)

write.csv(mesh$nodes, file.path(directory,"mesh_vertices.csv"))
write.csv(mesh$triangles, file.path(directory, "mesh_elements.csv"))
write.csv(mesh$nodesmarker, file.path(directory, "mesh_boundary.csv"))
write.csv(data, file.path(directory, "sample.csv"))

# Proposals for the smoothing parameter
lambda <- 10^seq(from = 0, to = -4, by = -1)

# Compute the log-density estimate
solution <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda,
                   nfolds = 10, preprocess_method = "RightCV",
                   tol1 = 1e-4, nsimulations = 1000,
                   step_method = "Wolfe_Method", direction_method = "BFGS")

best_lambda_id <- match(de$lambda, lambda_proposal)

write.csv(de$f_init[,best_lambda_id], file.path(directory, "f_init.csv"))

# Define the spatial field by FEM expansion
FEMfunction <- FEM(coeff = solution$g, FEMbasis = FEMbasis)

# Create a finer regular mesh of the surface of the unit sphere
mesh.eval <- refine.by.splitting.mesh.2.5D(mesh = mesh)

# Set up the finite element basis
FEMbasis.eval <- create.FEM.basis(mesh = mesh.eval)

# True density
true_density <- dens.func.4(data = mesh.eval$nodes, mesh = mesh)  # function in helper_functions_data.R

write.csv(true_density, file.path(directory, "true_density.csv"))

# Evaluate the log-density estimate on the nodes of the finer mesh
log_density_est <- read.csv("outputs/lbfgs30_curved_log_density.csv")
FEMresult <- FEM(coeff = log_density_est$V0, FEMbasis = FEMbasis)
evaluation <- eval.FEM(FEM = FEMresult, locations = mesh.eval$nodes)

# Compute the density estimate on the nodes of the finer mesh
estimated_density <- exp(evaluation)

# True density
true_density <- dens.func.4(data = mesh.eval$nodes, mesh = mesh)  # function in helper_functions_data.R

m = min(true_density, estimated_density, na.rm = TRUE)
M = max(true_density, estimated_density, na.rm = TRUE)

# Plot
plot.density.2.5D(coeff = true_density, mesh = mesh.eval, min_range = m, max_range = M,
                  colorscale = "viridis") 

plot.density.2.5D(coeff = estimated_density, mesh = mesh.eval, min_range = m, max_range = M,
                  colorscale = "viridis")
