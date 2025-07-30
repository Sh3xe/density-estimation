library(fdaPDE)
library(viridis)
rm(list = ls())

source("scripts/helper_functions_data.R")
source("scripts/helper_functions_plot.R")

TOL <- 1e-5
STEP <- 1e-2
MAX_ITER <- 500
CV_K <- 5

file_name <- function(test_directory, optimizer) {
	return( paste( c("outputs/r", optimizer, test_directory, "log_density.csv"), collapse = "_") )
}

load_scenario_2d <- function(dir_path) {
	dir <- file.path("./data/", dir_path)

	vertices <- read.csv(file.path(dir, "mesh_vertices.csv"))[, c("V1", "V2")]
	elements <- read.csv(file.path(dir, "mesh_elements.csv"))[, c("V1", "V2", "V3")]
	boundary <- read.csv(file.path(dir, "mesh_boundary.csv"))[, c("V1")]
	sample <- read.csv(file.path(dir, "sample.csv"))[, c("V1", "V2")]

	mesh <-create.mesh.2D(nodes=vertices, nodesattributes = boundary, triangles=elements)

	return( list(mesh=mesh, sample=sample, dir=dir_path) )
}

test_scenario <- function(scenario, lambda_proposal, direction_method) {
	FEMbasis <- create.FEM.basis(mesh = scenario$mesh)

	start_time <- proc.time()
	de <- DE.FEM(
		data = scenario$sample,
		FEMbasis = FEMbasis,
		lambda <- lambda_proposal,
		nfolds = CV_K,
		tol1 = TOL,
		nsimulations = MAX_ITER,
		heatIter = 1,
		preprocess_method = "RightCV",
		stepProposals = STEP,
		step_method = "Wolfe_Method",
		direction_method = direction_method
	)

	elapsed_time <- proc.time() - start_time
	best_lambda_id <- match(de$lambda, lambda_proposal)

	return(list(
		elapsed_time=elapsed_time[1],
		cv_err=de$CV_err[best_lambda_id],
		lambda=de$lambda,
		g=de$g
	))
}

benchmark_scenario <- function(scenario, lambdas) {
	methods <- c("L-BFGS10", "Gradient", "ConjugateGradientFR", "ConjugateGradientPRP")

	df <- data.frame(row.names=c("method", "cv_err", "lambda", "duration"))

	for(method in methods) {
		res <- test_scenario(scenario, lambdas, method)
		write.csv(res$g, file=file_name(scenario$dir, method))
		df <- rbind(df, list( "method" = method, "cv_err" = res$cv_err, "lambda" = res$lambda, "duration" = res$elapsed_time ))
	}

	return( df )
}

# Benchmark DE
lambdas_gs <- 10^seq(from = -1, to = -5, by = -0.5)
gs <- benchmark_scenario(load_scenario_2d("gaussian_square"), lambdas_gs)
write.csv(gs, "outputs/r_bmrk_gaussian_square.csv")

lambdas_is <- 10^seq(from = -1, to = -5, by = -0.5)
is <- benchmark_scenario(load_scenario_2d("infections_southampton"), lambdas_is)
write.csv(is, "outputs/r_bmrk_infection_southampton.csv")

lambdas_hs <- 10^seq(from = -1, to = -5, by = -1)
hs <- benchmark_scenario(load_scenario_2d("horseshoe"), lambdas_hs)
write.csv(hs, "outputs/r_bmrk_horseshoe.csv")