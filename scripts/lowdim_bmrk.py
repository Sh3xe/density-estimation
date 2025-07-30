import numpy as np
import pandas as pd
from time import time_ns
import scipy.optimize as opt
import pygad

MAX_ITER = 500
POP_SIZE = 100
TOL = 1e-5
STEP = 1e-2

DEFAULT_GAD_OPTS = {
	"num_generations": 50,
  "num_parents_mating": 5,   
	"random_seed": 0,
  "stop_criteria": "saturate_5",
	# "parent_selection_type": #"rank", "rws", "tournament"
	# "K_tournament": 2 # if parent_selection_type = tournament
  "keep_parents": -1, # could be 0 (keep no parent)
	"mutation_probability": 0.3,
	"mutation_type": "random",
  # keep_elitism=1
  # crossover_type="single_point"
  # crossover_probability=None
  # mutation_type="random"
  # mutation_by_replacement=False
	# mutation_percent_genes="default"
}

def duration_to_str(duration: int) -> str:
	if duration < 1e3:
		return f"{int(duration)}ns"
	elif duration < 1e6:
		return f"{int(duration * 1e-3)}\\textmu s"
	elif duration < 1e9:
		return f"{int(duration * 1e-6)}ms"
	elif duration < 1e12:
		return f"{int(duration * 1e-9)}s"
	else:
		return ">1000s"

def schaffer_f6(x: np.ndarray):
	r = x.dot(x)
	a = np.sin(np.sqrt(r))
	b = 1.0 + 0.001 * r
	return 0.5 + (a**2 - 0.5) / (b**2)

def rosenbrock(x: np.ndarray):
	a = 1.0
	b = 100.0
	return (a - x[0])**2 + b * (x[1] - x[0]**2)**2

def sphere(x: np.ndarray):
	return np.dot(x, x)

def schwefel(x: np.ndarray):
	s = 0.0
	for xi in x:
		if xi > 500 or xi < -500:
			s += 0.02 * xi * xi
		else:
			s += -xi * np.sin(np.sqrt(abs(xi)))
	return 418.9829 * x.size + s

def rastrigin(x: np.ndarray):
	s = 0.0
	for xi in x:
		if xi > 5.12 or xi < -5.12:
			s += 10.0 * xi**2
		else:
			s += xi**2 - 10.0 * np.cos(2 * np.pi * xi)
	return x.size * 10 + s

def benchmark_problem_scipy(fun, initial_points, method, method_min, title):
	if method == "L-BFGS-B":
		options = {"maxiter": MAX_ITER, "maxcor": 30}
	else:
		options = {"maxiter": MAX_ITER}

	# optimize
	results = pd.DataFrame(columns=["nit", "x_diff", "f_diff", "p_id", "duration_microsec"])
	for i in range(len(initial_points)):
		time_begin = time_ns()
		res = opt.minimize(
			fun=fun,
			x0=initial_points[i],
			method=method,
			tol=TOL,
			options=options
		)
		duration = (time_ns() - time_begin)

		results.loc[len(results)] = {
			"nit": res.nit,
			"x_diff": np.linalg.norm(res.x - method_min),
			"f_diff": np.linalg.norm(res.fun),
			"p_id": i,
			"duration_microsec": duration
		}	

	# Save benchmark as csv
	results.to_csv(f"outputs/py_bmrk_{method}_{title}_{len(initial_points)}.csv")
	# Print the result to a markdown table
	nit, nit_std = results["nit"].mean(), results["nit"].std()
	x_diff, x_diff_std = results["x_diff"].mean(), results["x_diff"].std()
	f_diff, f_diff_std = results["f_diff"].mean(), results["f_diff"].std()
	duration, duration_std = results["duration_microsec"].mean(), results["duration_microsec"].std()
	print(f"{title}|{nit:.4}|{nit_std:.4}|{x_diff:.2e}|{x_diff_std:.2e}|{f_diff:.2e}|{f_diff_std:.2e}|{duration_to_str(duration)}|{duration_to_str(duration_std)}")

def benchmark_problem_gad(fun, initial_points, method_opts, method_min, title):
	# optimize
	results = pd.DataFrame(columns=["nit", "x_diff", "f_diff", "p_id", "duration_microsec"])
	fitness_func = lambda ga, sol, idx: fun(sol)
	mutation_bound = sum([np.linalg.norm(v) for v in initial_points]) / len(initial_points)
	for i in range(len(initial_points)):
		initial_pop = [np.copy(initial_points[i]) for _ in range(POP_SIZE)]
		time_begin = time_ns()
		ga_instance = pygad.GA(
			fitness_func=fitness_func,
			initial_population=initial_pop,
			random_mutation_min_val = -mutation_bound,
			random_mutation_max_val = mutation_bound,
			**method_opts	)
		ga_instance.run()
		duration = (time_ns() - time_begin)
		solution, solution_fitness, _ = ga_instance.best_solution()

		results.loc[len(results)] = {
			"nit": ga_instance.best_solution_generation,
			"x_diff": np.linalg.norm(solution - method_min),
			"f_diff": np.linalg.norm(solution_fitness),
			"p_id": i,
			"duration_microsec": duration
		}	

	# Save benchmark as csv
	results.to_csv(f"outputs/py_bmrk_gad_{title}_{len(initial_points)}.csv")
	# Print the result to a markdown table
	nit, nit_std = results["nit"].mean(), results["nit"].std()
	x_diff, x_diff_std = results["x_diff"].mean(), results["x_diff"].std()
	f_diff, f_diff_std = results["f_diff"].mean(), results["f_diff"].std()
	duration, duration_std = results["duration_microsec"].mean(), results["duration_microsec"].std()
	print(f"{title}|{nit:.4}|{nit_std:.4}|{x_diff:.2e}|{x_diff_std:.2e}|{f_diff:.2e}|{f_diff_std:.2e}|{duration_to_str(duration)}|{duration_to_str(duration_std)}")

def test_method(scipy_method, gad_dict):
	print(f"### `{scipy_method}`\n")
	print("| Method | n_iter_mean | n_iter_std | x_diff_mean | x_diff_std | f_diff_mean | f_diff_std | duration_mean | duration_std ")
	print("|-|-|-|-|-|-|-|-|-|")

	init_sphere_2d = np.loadtxt("data/lowdim_inits/2_sphere.csv", delimiter=",")
	init_sphere_10d = np.loadtxt("data/lowdim_inits/10_sphere.csv", delimiter=",")
	init_sphere_30d = np.loadtxt("data/lowdim_inits/30_sphere.csv", delimiter=",")
	init_schwefel_2d = np.loadtxt("data/lowdim_inits/2_schwefel.csv", delimiter=",")
	init_schwefel_10d = np.loadtxt("data/lowdim_inits/10_schwefel.csv", delimiter=",")
	init_schwefel_30d = np.loadtxt("data/lowdim_inits/30_schwefel.csv", delimiter=",")
	init_rastrigin_2d = np.loadtxt("data/lowdim_inits/2_rastrigin.csv", delimiter=",")
	init_rastrigin_10d = np.loadtxt("data/lowdim_inits/10_rastrigin.csv", delimiter=",")
	init_rastrigin_30d = np.loadtxt("data/lowdim_inits/30_rastrigin.csv", delimiter=",")
	init_schaffer_f6 = np.loadtxt("data/lowdim_inits/schaffer_f6.csv", delimiter=",")
	init_rosenbrock = np.loadtxt("data/lowdim_inits/rosenbrock.csv", delimiter=",")

	if scipy_method != None:
		benchmark_problem_scipy(sphere, init_sphere_2d, scipy_method, np.zeros(2), "sphere_2d")
		benchmark_problem_scipy(sphere, init_sphere_10d, scipy_method, np.zeros(10), "sphere_10d")
		benchmark_problem_scipy(sphere, init_sphere_30d, scipy_method, np.zeros(30), "sphere_30d")
		benchmark_problem_scipy(schwefel, init_schwefel_2d, scipy_method, np.ones(2) * 420.9687, "schwefel_2d")
		benchmark_problem_scipy(schwefel, init_schwefel_10d, scipy_method, np.ones(10) * 420.9687, "schwefel_10d")
		benchmark_problem_scipy(schwefel, init_schwefel_30d, scipy_method, np.ones(30) * 420.9687, "schwefel_30d")
		benchmark_problem_scipy(rastrigin, init_rastrigin_2d, scipy_method, np.zeros(2), "rastrigin_2d")
		benchmark_problem_scipy(rastrigin, init_rastrigin_10d, scipy_method, np.zeros(10), "rastrigin_10d")
		benchmark_problem_scipy(rastrigin, init_rastrigin_30d, scipy_method, np.zeros(30),"rastrigin_30d")
		benchmark_problem_scipy(schaffer_f6, init_schaffer_f6, scipy_method, np.zeros(2), "schaffer_f6")
		benchmark_problem_scipy(rosenbrock, init_rosenbrock, scipy_method, np.ones(2), "rosenbrock")
	elif gad_dict != None:
		benchmark_problem_gad(sphere, init_sphere_2d, gad_dict, np.zeros(2), "Sphere 2D")
		benchmark_problem_gad(sphere, init_sphere_10d, gad_dict, np.zeros(10), "Sphere 10D")
		benchmark_problem_gad(sphere, init_sphere_30d, gad_dict, np.zeros(30), "Sphere 30D")
		benchmark_problem_gad(schwefel, init_schwefel_2d, gad_dict, np.ones(2) * 420.9687, "schwefel_2d")
		benchmark_problem_gad(schwefel, init_schwefel_10d, gad_dict, np.ones(10) * 420.9687, "schwefel_10d")
		benchmark_problem_gad(schwefel, init_schwefel_30d, gad_dict, np.ones(30) * 420.9687, "schwefel_30d")
		benchmark_problem_gad(rastrigin, init_rastrigin_2d, gad_dict, np.zeros(2), "rastrigin_2d")
		benchmark_problem_gad(rastrigin, init_rastrigin_10d, gad_dict, np.zeros(10), "rastrigin_10d")
		benchmark_problem_gad(rastrigin, init_rastrigin_30d, gad_dict, np.zeros(30),"rastrigin_30d")
		benchmark_problem_gad(schaffer_f6, init_schaffer_f6, gad_dict, np.zeros(2), "schaffer_f6")
		benchmark_problem_gad(rosenbrock, init_rosenbrock, gad_dict, np.ones(2), "rosenbrock")

if __name__ == "__main__":
	test_method("L-BFGS-B", None)
	test_method("Nelder-Mead", None)
	test_method("CG", None)