library(fdaPDE)
library(viridis)
rm(list = ls())

# This script is located under [repo/scripts] and the data is located at [repo/data] 
setwd(file.path(getwd(), ".."))

# Import helper functions to sample and visualize data
source("scripts/helper_functions_plot.R")
source("scripts/helper_functions_data.R")

data_dir <- "data/gaussian_2d"
nodes <- read.csv(file.path(data_dir,"mesh_vertices.csv"))[,c(2,3)]
triangles <- read.csv(file.path(data_dir, "mesh_elements.csv"))[,c(2,3,4)]
nodesmarker <- read.csv(file.path(data_dir, "mesh_boundary.csv"))[,2]

# Create a regular mesh of the unit square
mesh <- create.mesh.2D(nodes = nodes, triangles = triangles, nodesattributes = nodesmarker)
# Set up the finite element basis
FEMbasis <- create.FEM.basis(mesh = mesh)

# Load the data
output_dir <- "data/gaussian_2d"
f_init <- read.csv(file.path(output_dir, "f_init.csv"))
sample <- read.csv(file.path(output_dir, "sample.csv"))
log_dens <- read.csv(file.path(output_dir, "log_dens.csv"))

# Plot
par(mfrow = c(1,2), mai = c(0.5,0.25,0.5,0.5))
x_plot <- seq(from=0.0, to=1.0, length.out = 100)
y_plot <- seq(from=0.0, to=1.0, length.out = 100)
grid_plot <- expand.grid(x_plot, y_plot)
fem_log_dens <- FEM(coeff = log_dens$x, FEMbasis = FEMbasis)
log_dens_eval <- eval.FEM(FEM = fem_log_dens, grid_plot)
dens_eval <- exp(log_dens_eval)

plot.density.2D(
	X = x_plot, Y = y_plot, Z = dens_eval,
	colorscale = "viridis",
	boundary = matrix(
		data = c(0,0,1,0,1,1,0,1,0,0),
		ncol = 2,
		byrow = TRUE
	)
)

# Generate a high res mesh
hr_grid <- expand.grid(x_plot, y_plot)
bounds <- hr_grid[(hr_grid$Var1 %in% c(0, 1)) | (hr_grid$Var2 %in% c(0, 1)), ]
hr_mesh <- create.mesh.2D(nodes=hr_grid)

plot.density.2D(
	X = x_plot, Y = y_plot, Z = dens.func.1(data=hr_mesh$nodes),
	colorscale = "viridis",
	boundary = matrix(
		data = c(0,0,1,0,1,1,0,1,0,0),
		ncol = 2,
		byrow = TRUE
	)
)
