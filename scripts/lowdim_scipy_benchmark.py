import numpy as np
import pandas as pd
from time import time_ns
import scipy.optimize as opt

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

def benchmark_problem(fun, initial_points, method, method_min, title):
	if method == "L-BFGS-B":
		options = {"maxiter": 500, "maxcor": 30}
	else:
		options = {"maxiter": 500}

	# optimize
	time_begin = time_ns()
	results = pd.DataFrame(columns=["nit", "x_diff", "f_diff", "p_id"])
	for i in range(len(initial_points)):
		res = opt.minimize(
			fun=fun,
			x0=initial_points[i],
			method=method,
			tol=1e-5,
			options=options
		)

		results.loc[len(results)] = {
			"nit": res.nit,
			"x_diff": np.linalg.norm(res.x - method_min),
			"f_diff": np.linalg.norm(res.fun),
			"p_id": i
		}

	duration = (time_ns() - time_begin) / len(initial_points)

	nit = results["nit"].mean()
	x_diff = results["x_diff"].mean()
	x_diff_std = results["x_diff"].std()
	f_diff = results["f_diff"].mean()
	f_diff_std = results["f_diff"].std()
	print(f"{title}| {nit:.4} | {x_diff:.2e} | {x_diff_std:.2e} | {f_diff:.2e} |  {f_diff_std:.2e}| {duration_to_str(duration)}")

def test_method(method):

	print(f"### `{method}`\n")
	print("|Method|Iterations| x_diff | x_diff_std | f_diff | f_diff_std | time |")
	print("|-|-|-|-|-|-|-|")

	init_sphere_2d = np.loadtxt("data/lowdim_inits/2_sphere.csv", delimiter=",")
	benchmark_problem(sphere, init_sphere_2d, method, np.zeros(2), "Sphere 2D")
	init_sphere_10d = np.loadtxt("data/lowdim_inits/10_sphere.csv", delimiter=",")
	benchmark_problem(sphere, init_sphere_10d, method, np.zeros(10), "Sphere 10D")
	init_sphere_30d = np.loadtxt("data/lowdim_inits/30_sphere.csv", delimiter=",")
	benchmark_problem(sphere, init_sphere_30d, method, np.zeros(30), "Sphere 30D")

	init_schwefel_2d = np.loadtxt("data/lowdim_inits/2_schwefel.csv", delimiter=",")
	benchmark_problem(schwefel, init_schwefel_2d, method, np.ones(2) * 420.9687, "Schwefel 2D")
	init_schwefel_10d = np.loadtxt("data/lowdim_inits/10_schwefel.csv", delimiter=",")
	benchmark_problem(schwefel, init_schwefel_10d, method, np.ones(10) * 420.9687, "Schwefel 10D")
	init_schwefel_30d = np.loadtxt("data/lowdim_inits/30_schwefel.csv", delimiter=",")
	benchmark_problem(schwefel, init_schwefel_30d, method, np.ones(30) * 420.9687, "Schwefel 30D")


	init_rastrigin_2d = np.loadtxt("data/lowdim_inits/2_rastrigin.csv", delimiter=",")
	benchmark_problem(rastrigin, init_rastrigin_2d, method, np.zeros(2), "Rastrigin 2D")
	init_rastrigin_10d = np.loadtxt("data/lowdim_inits/10_rastrigin.csv", delimiter=",")
	benchmark_problem(rastrigin, init_rastrigin_10d, method, np.zeros(10), "Rastrigin 10D")
	init_rastrigin_30d = np.loadtxt("data/lowdim_inits/30_rastrigin.csv", delimiter=",")
	benchmark_problem(rastrigin, init_rastrigin_30d, method, np.zeros(30),"Rastrigin 30D")

	init_schaffer_f6 = np.loadtxt("data/lowdim_inits/schaffer_f6.csv", delimiter=",")
	benchmark_problem(schaffer_f6, init_schaffer_f6, method, np.zeros(2), "Schaffer F6")

	init_rosenbrock = np.loadtxt("data/lowdim_inits/rosenbrock.csv", delimiter=",")
	benchmark_problem(rosenbrock, init_rosenbrock, method, np.ones(2), "Rosenbrock")

if __name__ == "__main__":
	test_method("L-BFGS-B")
	test_method("Nelder-Mead")
	test_method("CG")