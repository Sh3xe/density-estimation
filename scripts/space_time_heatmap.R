library(fdaPDE)

dir <- "data/infections_southampton"
mesh_vertices <- read.csv(file.path(dir, "mesh_vertices.csv"))[,c("V1", "V2")]
mesh_boundaries <- read.csv(file.path(dir, "mesh_boundary.csv"))$V1
mesh_elements <- read.csv(file.path(dir, "mesh_elements.csv"))[,c("V1", "V2", "V3")]

sample_space <- read.csv(file.path(dir, "sample.csv"))[,c("V1", "V2")]
sample_time <- read.csv(file.path(dir, "sample_time.csv"))$x

mesh <- create.mesh.2D(nodes = mesh_vertices, nodesattributes = mesh_boundaries, triangles = mesh_elements)
fem_basis <- create.FEM.basis(mesh = mesh)

mesh_time <- seq(from=min(sample_time), to=max(sample_time), length.out = 9)

CV_K <- 5
MAX_ITER <- 200
TOL <- 1e-3

space_values <- 10^seq(from=-2, to=2, by=0.25)
time_values <- 10^seq(from=-4, to=1, by=0.25)

de_res <- DE.FEM.time(
	sample_space, sample_time,
	fem_basis, mesh_time,
	space_values, time_values,
	step_method = "Wolfe_Method",
	direction_method = "L-BFGS10",
	preprocess_method = "RightCV",
	nfolds = CV_K,
	nsimulations = MAX_ITER,
	tol1 = TOL, tol2 = TOL
)

write.csv(space_values, "outputs/i_st_hm_space.csv")
write.csv(time_values, "outputs/i_st_hm_time.csv")
write.csv(de_res$CV_err, "outputs/i_st_hm_cverr.csv")