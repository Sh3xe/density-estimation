import numpy as np
from time import time_ns
import scipy.optimize as opt

def duration_to_str(duration: int) -> str:
	if duration < 1e3:
		return f"{int(duration)}ns"
	elif duration < 1e6:
		return f"{int(duration * 1e-3)}μs"
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

def benchmark_problem(fun, x0, method, method_min, title):
	if method == "L-BFGS-B":
		options = {"maxiter": 500, "maxcor": 30}
	else:
		options = {"maxiter": 500}

	optimize_fun = lambda: opt.minimize(
		fun=fun,
		x0=x0,
		method=method,
		tol=1e-5,
		options=options
	)
	
	time_begin = time_ns()
	for _ in range(100):
		optimize_fun()
	
	duration = (time_ns() - time_begin) / 100
	
	res = optimize_fun()

	x_diff = np.linalg.norm(res.x - method_min)
	print(f"{title}| {res.nit} | {x_diff:.2e} | {np.linalg.norm(res.fun):.2e} | {duration_to_str(duration)}")

def test_method(method):

	print(f"### `{method}`\n")
	print("|Method|Iterations| x-x* l2 err | f-f* l2 err | time |")
	print("|-|-|-|-|-|")
	benchmark_problem(sphere, np.ones(2), method, np.zeros(2), "Sphere 2D")
	benchmark_problem(sphere, np.ones(30), method, np.zeros(30), "Sphere 30D")
	benchmark_problem(schwefel, np.ones(2) * 20, method, np.ones(2) * 420.9687, "Schwefel 2D")
	benchmark_problem(schwefel, np.ones(10) * 20, method, np.ones(10) * 420.9687, "Schwefel 10D")
	benchmark_problem(rastrigin, np.ones(2) * 3, method, np.zeros(2), "Rastrigin 2D")
	benchmark_problem(rastrigin, np.ones(30) * 3, method, np.zeros(30),"Rastrigin 30D")
	benchmark_problem(schaffer_f6, np.ones(2) * 5, method, np.zeros(2), "Schaffer F6")
	benchmark_problem(rosenbrock, np.array([-1.2, 1]), method, np.ones(2), "Rosenbrock")

if __name__ == "__main__":
	test_method("L-BFGS-B")
	test_method("Nelder-Mead")
	test_method("CG")