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

split_dataset <- function(dataset, k, i) {
  slice_size <- nrow(dataset) / k
  slice_index_begin <- min(i * slice_size, nrow(dataset) - slice_size)

  testing <- dataset[(slice_index_begin + 1):(slice_index_begin + slice_size), ]

  training_1 <- dataset[1:slice_index_begin, ]
	if ((slice_index_begin + slice_size) < nrow(dataset)) {
    training_2 <- dataset[(slice_index_begin + slice_size + 1):nrow(dataset), ]
  } else {
    training_2 <- matrix(nrow = 0, ncol = ncol(dataset))  # Empty matrix with the same number of columns
  }


  training <- rbind(training_1, training_2)

  return(list(training = training, testing = testing))
}

load_scenario_2d <- function(dir_path, scenario_title) {
	dir <- file.path("./data/", dir_path)

	vertices <- read.csv(file.path(dir, "mesh_vertices.csv"))[, c("V1", "V2")]
	elements <- read.csv(file.path(dir, "mesh_elements.csv"))[, c("V1", "V2", "V3")]
	boundary <- read.csv(file.path(dir, "mesh_boundary.csv"))[, c("V1")]
	sample <- read.csv(file.path(dir, "sample.csv"))[, c("V1", "V2")]

	mesh <-create.mesh.2D(nodes=vertices, nodesattributes = boundary, triangles=elements, order=1)

	return( list(mesh=mesh, sample=sample, dir=dir_path, title=scenario_title) )
}

test_scenario <- function(scenario, lambda_proposal, direction_method) {
	FEMbasis <- create.FEM.basis(mesh = scenario$mesh)

	start_time <- proc.time()
	de <- DE.FEM(
		data = scenario$sample,
		FEMbasis = FEMbasis,
		lambda <- lambda_proposal,
		# nfolds = CV_K,
		tol1 = TOL,
		nsimulations = MAX_ITER,
		# preprocess_method = "RightCV",
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
	methods <- c("L-BFGS10", "BFGS", "ConjugateGradientPRP")

	df <- data.frame(row.names=c("method", "cv_err", "lambda", "duration"))

	for(method in methods) {
		res <- test_scenario(scenario, lambdas, method)
		write.csv(res$g, file=file_name(scenario$dir, method))
		df <- rbind(df, list( "method" = method, "cv_err" = res$cv_err, "lambda" = res$lambda, "duration" = res$elapsed_time ))
	}

	return( df )
}

cv_error <- function(log_density, scenario, cv_k) {
	FEMbasis <- create.FEM.basis(mesh = scenario$mesh)

	err <- 0.0
	for (k in 1:cv_k) {
		res <- split_dataset(scenario$sample, cv_k, k)
		fem_log_dens <- FEM(coeff = log_density, FEMbasis = FEMbasis)
		log_dens_eval <- eval.FEM(FEM = fem_log_dens, res$testing)
		err <- err - (sum(log_dens_eval) / length(log_dens_eval))
	}

	return (err / cv_k)
}

compute_cv_errors <- function(language_str, opt_method, scenario, lambda) {
	files <- list.files("cv_csv")
	n_errs <- 0
	
	for(filename in files) {
		parts <- strsplit(strsplit(filename, ".csv")[[1]], "#")[[1]]

		file_language_str <- parts[1]
		file_i_str <- parts[2]
		file_opt_method <- parts[3]
		file_lambda_str <- parts[4]
		file_scenario_title <- parts[5]

		i <- as.integer(file_i_str)
		file_lambda <- as.double(file_lambda_str)
		if( (file_language_str == language_str) && (file_opt_method == opt_method) && (file_scenario_title == scenario$title) && (abs(file_lambda - lambda) < 1e-5) ) {
			log_dens <- read.csv(file.path("cv_csv", filename))
			if(file_language_str == "r") {
				log_dens_vec <- log_dens$x
			} else {
				log_dens_vec <- log_dens$V0
			}
			cv_err <- cv_error(log_dens_vec, scenario, 5)
			err_tot <- cv_err
			n_errs <- n_errs + 1
		}
	}

	return (err_tot / n_errs)
}

comp_densities <- function(scenario, r_title, cpp_title) {
	FEMbasis <- create.FEM.basis(mesh = scenario$mesh)
	r_log_dens <- read.csv(r_title)$x
	cpp_log_dens <- read.csv(cpp_title)$V0

	l2_norm <- r_log_dens - cpp_log_dens
	fem_diff <- FEM(coeff = l2_norm, FEMbasis = FEMbasis)
	
	plot(fem_diff)
}

output_scenario_cv_densities <- function(scenario, lambda_proposal, direction_method, cv_k) {
	FEMbasis <- create.FEM.basis(mesh = scenario$mesh)

	for(cv_i in 0:(cv_k-1)) {
		dataset <- split_dataset(scenario$sample, cv_k, cv_i)

		for(lambda in lambda_proposal) {
			de <- DE.FEM(
				data = dataset$training,
				FEMbasis = FEMbasis,
				lambda <- lambda,
				tol1 = TOL,
				nsimulations = MAX_ITER,
				stepProposals = STEP,
				step_method = "Wolfe_Method",
				direction_method = direction_method
			)

			file_name <- paste0("cv_csv/r#", cv_i, "#", opt_name, "#", lambda, "#", scenario$title, ".csv")
			write.csv(de$g, file_name)
		}
	}
}

compute_all_cv_errors <- function(cpp_opt_names, r_opt_names, scenarios, lambda_proposal) {
	results_df <- data.frame(
		r_opt_name = character(),
		cpp_opt_name = character(),
		lambda = numeric(),
		scenario = character(),
		r_cv_err = numeric(),
		cpp_cv_err = numeric(),
		cv_err_diff = numeric(),
		stringsAsFactors = FALSE
	)

	for(i_opt in 1:length(cpp_opt_names)) {
		r_opt_name <- r_opt_names[i_opt]
		cpp_opt_name <- cpp_opt_names[i_opt]
		for(lambda in lambda_proposal) {
			for(scenario in scenarios) {
				r_cv_err <- compute_cv_errors("r", r_opt_name, scenario, lambda)
				cpp_cv_err <- compute_cv_errors("cpp", cpp_opt_name, scenario, lambda)
				print(paste("r_cv_err:", r_cv_err, ", cpp_cv_err:", cpp_cv_err))
				results_df <- rbind(results_df, data.frame(
					r_opt_name = r_opt_name,
					cpp_opt_name = cpp_opt_name,
					lambda = lambda,
					scenario = scenario$title,
					r_cv_err = r_cv_err,
					cpp_cv_err = cpp_cv_err,
					cv_err_diff = abs(r_cv_err - cpp_cv_err)
				))
			}
		}
	}	
	write.csv(results_df, "cross_validation_errors.csv", row.names = FALSE)
}

compute_best_lambdas <- function(csv_table_path, cpp_opt_names, scenarios) {
  df <- read.csv(csv_table_path)

  for (scenario in scenarios) {
    for (cpp_opt_name in cpp_opt_names) {
      rows <- df[df$cpp_opt_name == cpp_opt_name & df$scenario == scenario$title, ]

      if (nrow(rows) > 0) {
        best_r <- rows[which.min(rows$r_cv_err), ]
				best_cpp <- rows[which.min(rows$cpp_cv_err), ]

				print(paste("scenario=", scenario$title, ", opt=",cpp_opt_name ,"best_r=", best_r$lambda, ", best_cpp=", best_cpp$lambda))
      } else {
        message(paste("No matching rows found for scenario:", scenario$title, "and cpp_opt_name:", cpp_opt_name))
      }
    }
  }
}


gaussian_square <- load_scenario_2d("gaussian_square", "gaussian_square")
horseshoe <- load_scenario_2d("horseshoe", "horseshoe")
infections_southampton <- load_scenario_2d("infections_southampton", "infections_southampton")
lambda_proposal <- 10^seq(from = 2, to = -2, by = -1.0)
r_opt_names <- c("BFGS", "L-BFGS10", "ConjugateGradientPRP")
cpp_opt_names <- c("BFGS", "LBFGS10", "cg_pr")
scenarios <- list(gaussian_square, horseshoe, infections_southampton)
cv_k <- 5

compute_best_lambdas("cross_validation_errors.csv", cpp_opt_names, scenarios)

# compute_all_cv_errors(cpp_opt_names, r_opt_names, scenarios, lambda_proposal)

# for(opt_name in r_opt_names) {
	# output_scenario_cv_densities(gaussian_square, lambda_proposal, opt_name, cv_k)
	# output_scenario_cv_densities(horseshoe, lambda_proposal, opt_name, cv_k)
	# output_scenario_cv_densities(infections_southampton, lambda_proposal, opt_name, cv_k)
# }

# comp_densities(infections_southampton, "outputs/r_ConjugateGradientPRP_infections_southampton_log_density.csv", "outputs/cpp_cg_pr_infections_southampton_log_density.csv")
# compute_cv_errors("cpp", "LBFGS10", load_scenario_2d("gaussian_square", "gaussian_square"))
# compute_cv_errors("cpp", "BFGS", load_scenario_2d("gaussian_square", "gaussian_square"))
# compute_cv_errors("cpp", "cg_pr", load_scenario_2d("infections_southampton", "infections_southampton"))

# Benchmark DE
# gs <- benchmark_scenario(load_scenario_2d("gaussian_square", "gaussian_square"), c(0.1))
# write.csv(gs, "outputs/r_bmrk_gaussian_square.csv")

# is <- benchmark_scenario(load_scenario_2d("infections_southampton", "infections_southampton"), c(0.1))
# write.csv(is, "outputs/r_bmrk_infection_southampton.csv")

# hs <- benchmark_scenario(load_scenario_2d("horseshoe", "horseshoe"), c(0.1))
# write.csv(hs, "outputs/r_bmrk_horseshoe.csv")

