import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import json
import imageio

def schaffer_f6(x, y):
	r = x**2 + y**2
	a = np.sin(np.sqrt(r))
	b = 1.0 + 0.001 * r
	return 0.5 + (a**2 - 0.5) / (b**2)

def rosenbrock(x, y):
	a = 1.0
	b = 100.0
	return (a - x)**2 + b * (y - x**2)**2

def sphere(x, y):
	return x**2 + y**2

def schwefel(x, y):
	s = 0.0
	for xi in [x, y]:
		if xi > 500 or xi < -500:
			s += 0.02 * xi * xi
		else:
			s += -xi * np.sin(np.sqrt(abs(xi)))
	return 418.9829 * 2 + s

def rastrigin(x, y):
	s = 0.0
	for xi in [x, y]:
		if xi > 5.12 or xi < -5.12:
			s += 10.0 * xi**2
		else:
			s += xi**2 - 10.0 * np.cos(2 * np.pi * xi)
	return 2 * 10 + s

def plot_function_graph():
	x = np.arange(-2, 2, 0.1)
	y = np.arange(-2, 2, 0.1)
	X, Y = np.meshgrid(x, y)
	Z = np.vectorize(rosenbrock)(X, Y)

	fig = plt.figure(figsize=(10, 7))
	ax = fig.add_subplot(111, projection='3d')
	surf = ax.plot_surface(X, Y, Z, cmap='RdBu', rstride=1,cstride=1, edgecolor='none', shade=True)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('rosenbrock')
	ax.set_proj_type('ortho')
	plt.show()

def box_plot_diff(py_outs, cpp_outs, method_title, title, ax, use_log=False):
	py_csv = list(filter(lambda x: method_title in x, py_outs))
	cpp_csv = list(filter(lambda x: method_title in x, cpp_outs))

	if len(py_csv) != 1 or len(cpp_csv) != 1:
		print("Cannot find csv files")
		return
	
	py_df = pd.read_csv("outputs/" + py_csv[0])
	cpp_df = pd.read_csv("outputs/" + cpp_csv[0])
	
	data = [cpp_df["x_diff"].to_numpy(), py_df["x_diff"].to_numpy()]
	ax.set_ylabel('x_diff')
	if use_log:
		ax.set_yscale("log")
	bplot = ax.boxplot(
		data,
		patch_artist=True,
		tick_labels=["fdaPDE", "scipy"]
	)

	bplot["boxes"][0].set_facecolor("#EF8A87")
	bplot["boxes"][1].set_facecolor("#85A2EB")

	ax.set_title(title)

def box_plot_lowdim(py_name, cpp_name, fig_title):
	files = os.listdir("outputs")
	outputs = list(filter(lambda x: ".csv" in x, files))
	py_outs = list(filter(lambda x: "py" in x and py_name in x, outputs))
	cpp_outs = list(filter(lambda x: "cpp" in x and cpp_name in x, outputs))

	test_cases = ["sphere", "rastrigin", "schwefel"]
	dimensions = [2, 10, 30]
	other_test_cases = ["rosenbrock", "schaffer_f6"]
	use_logs = [True, False, False]

	for j in range(len(test_cases)):
		test_case = test_cases[j]
		fig, axes = plt.subplots(1, 3, figsize=(18, 6))
		fig.suptitle(fig_title)
		for i, dim in enumerate(dimensions):
			method_title = f"{test_case}_{dim}d"
			title = f"{test_case.capitalize()} {dim}D"
			box_plot_diff(py_outs, cpp_outs, method_title, title, axes[i], use_logs[j])
		plt.tight_layout()
		plt.title(title)
		plt.savefig(f"figures/{test_case}_{cpp_name}_comparison.png")
		plt.close()

	for test_case in other_test_cases:
		fig, ax = plt.subplots(figsize=(6, 6))
		fig.suptitle(fig_title)
		box_plot_diff(py_outs, cpp_outs, test_case, test_case.capitalize(), ax)
		plt.yscale("linear")
		plt.tight_layout()
		plt.title(title)
		plt.savefig(f"figures/{test_case}_{cpp_name}_comparison.png")
		plt.close()

def load_population_data(file_path):
	with open(file_path, "r") as file:
		data = json.load(file)
	return data

def plot_population(population, time_frame):
	MAX_VAL = 500
	x = np.linspace(-MAX_VAL, MAX_VAL, MAX_VAL)
	y = np.linspace(-MAX_VAL, MAX_VAL, MAX_VAL)
	X, Y = np.meshgrid(x, y)
	Z = np.array([[schwefel(xi, yi) for xi in x] for yi in y])

	plt.imshow(Z, extent=[-MAX_VAL, MAX_VAL, -MAX_VAL, MAX_VAL], origin="lower", cmap="viridis", alpha=0.5)
	plt.colorbar(label="Schwefel Function Value")
	plt.scatter(*zip(*population), alpha=0.7, color="red")
	plt.title(f"Population at Time Frame {time_frame}")
	plt.xlabel("X Coordinate")
	plt.ylabel("Y Coordinate")
	plt.xlim(-MAX_VAL, MAX_VAL)
	plt.ylim(-MAX_VAL, MAX_VAL)
	plt.grid(True)

def create_population_plots(data):
	images = []
	for time_frame, population in enumerate(data):
		plt.figure()
		plot_population(population, time_frame + 1)
		plt.savefig(f"population_{time_frame}.png")
		images.append(imageio.imread(f"population_{time_frame}.png"))
	plt.close()
	return images

def create_gif(images, output_path):
	imageio.mimsave(output_path, images, fps=10)

def animation_genetic():
	file_path = "population.json"
	output_gif_path = "outputs/population_animation.gif"
	data = load_population_data(file_path)
	images = create_population_plots(data)
	create_gif(images, output_gif_path)

if __name__ == "__main__":
	box_plot_lowdim("L-BFGS-B", "lbfgs30", "L-BFGS-30")
	box_plot_lowdim("Nelder-Mead", "nelder_mead", "Nelder-Mead")
	box_plot_lowdim("CG", "cg_fr", "Conjugate Gradient")