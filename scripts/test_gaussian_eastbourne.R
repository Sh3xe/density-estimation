library(fdaPDE)
library(viridis)
rm(list = ls())
directory <- "data/gaussian_eastbourne"

# Import helper functions to sample and visualize data
source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")

# Define the domain
mesh = mesh.5()
mesh_linnet = as.linnet(mesh)
mesh$nodesmarkers[mesh$nodesmarkers == TRUE] <- 1
mesh$nodesmarkers[mesh$nodesmarkers == FALSE] <- 0

# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Generate the data
pp <- rpoislpp(lambda = intens.func.5, L = mesh_linnet, sigma = 0.2)
data <- cbind(pp$data$x, pp$data$y)

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

mesh.eval <- refine.mesh.1.5D(mesh = mesh, delta = 0.05)
mesh.eval_linnet <- as.linnet(mesh.eval)
FEMbasis.eval <- create.FEM.basis(mesh = mesh.eval)

# Calculate the true density function
true_density <- intens.func.5(
	x = mesh.eval$nodes[,1],
	y = mesh.eval$nodes[,2],
	L = mesh.eval_linnet, sigma = 0.2
) / integral.dens.func.5(f = FEMbasis)


write.csv(true_density, file.path(directory, "true_density.csv"))

plot.sample.1.5D(data = data, mesh = mesh, linewidth = 1.5)


# Plot
FEMfunction <- FEM(coeff = de$g, FEMbasis = FEMbasis)

mesh.eval <- refine.mesh.1.5D(mesh = mesh, delta = 0.05)
mesh.eval_linnet <- as.linnet(mesh.eval)
FEMbasis.eval <- create.FEM.basis(mesh = mesh.eval)

evaluation <- eval.FEM(FEM = FEMfunction, locations = mesh.eval$nodes)
estimated_density <- exp(evaluation)

# True density
true_density <- intens.func.5(x = mesh.eval$nodes[,1], y = mesh.eval$nodes[,2],
                              L = mesh.eval_linnet, sigma = 0.2) / integral.dens.func.5(f = FEMbasis)

# Plot
plot.density.1.5D(coeff = true_density, mesh = mesh.eval, colorscale = "jet.col",
                  linewidth = 1.5) +
plot.density.1.5D(coeff = estimated_density, mesh = mesh.eval, colorscale = "jet.col",
                  linewidth = 1.5)
