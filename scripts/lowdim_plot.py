import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import pandas as pd
import json
import imageio

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
	# cpp_df = cpp_df[cpp_df["x_diff"] < 1e3]
	
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
	py_outs = list(filter(lambda x: "py" in x and py_name + "_" in x, outputs))
	cpp_outs = list(filter(lambda x: "cpp" in x and cpp_name  + "_" in x, outputs))

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

	plt.imshow(Z, extent=[-MAX_VAL, MAX_VAL, -MAX_VAL, MAX_VAL], origin="lower",  alpha=0.5)
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

def plot_initial_points(function, bound, csv_name):
	# Find initial points
	df = pd.read_csv(f"data/lowdim_inits/{csv_name}.csv", names=["x", "y"])

	# Plots function
	x = np.linspace(-bound, bound, 100)
	y = np.linspace(-bound, bound, 100)
	X, Y = np.meshgrid(x, y)
	Z = np.array([[function(xi, yi) for xi in x] for yi in y])

	plt.imshow(Z, extent=[-bound, bound, -bound, bound], origin="lower", cmap="RdBu", alpha=0.5)
	# plt.colorbar(label="Schwefel Function Value")
	plt.scatter(x=df["x"], y=df["y"], alpha=0.7, color="#000000")
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.xlim(-bound, bound)
	plt.ylim(-bound, bound)
	plt.grid(True)
	plt.title(csv_name)
	plt.tight_layout()
	plt.savefig(f"figures/{csv_name}_init_point.png")
	plt.close()

def output_diff_table_line(function_name, py_outs, cpp_outs):
	py_csv = list(filter(lambda x: function_name in x, py_outs))
	cpp_csv = list(filter(lambda x: function_name in x, cpp_outs))

	if len(py_csv) != 1 or len(cpp_csv) != 1:
		print("Cannot find csv files")
		return
	
	py_df = pd.read_csv("outputs/" + py_csv[0])
	cpp_df = pd.read_csv("outputs/" + cpp_csv[0])

	s = ""

	s += "{} (fdaPDE)|{:.1f} \\textpm {:.1f} |{:.2e} \\textpm {:.2e}|{:.2e}\\textpm {:.2e}|{}\\textpm {}\n".format(
		function_name.replace("_", "\\_"),
		cpp_df["nit"].mean(), cpp_df["nit"].std(),
		cpp_df["x_diff"].mean(), cpp_df["x_diff"].std(),
		cpp_df["f_diff"].mean(), cpp_df["f_diff"].std(),
		duration_to_str(cpp_df["duration_microsec"].mean()), duration_to_str(cpp_df["duration_microsec"].std()),
		)
	
	s += "{} (scipy)|{:.1f} \\textpm {:.1f} |{:.2e} \\textpm {:.2e}|{:.2e}\\textpm {:.2e}|{}\\textpm {}\n".format(
		function_name.replace("_", "\\_"),
		py_df["nit"].mean(), py_df["nit"].std(),
		py_df["x_diff"].mean(), py_df["x_diff"].std(),
		py_df["f_diff"].mean(), py_df["f_diff"].std(),
		duration_to_str(py_df["duration_microsec"].mean()), duration_to_str(py_df["duration_microsec"].std()),
		)
	return s

def output_diff_table(py_name, cpp_name):
	files = os.listdir("outputs")
	outputs = list(filter(lambda x: ".csv" in x, files))
	py_outs = list(filter(lambda x: "py" in x and py_name + "_" in x, outputs))
	cpp_outs = list(filter(lambda x: "cpp" in x and cpp_name + "_" in x, outputs))

	s = "Method|Iters|x\\_diff|f\\_diff|duration\n"
	s += "-|-|-|-|-\n"

	s += output_diff_table_line("sphere_2d", py_outs, cpp_outs)
	s += output_diff_table_line("sphere_10d", py_outs, cpp_outs)
	s += output_diff_table_line("sphere_30d", py_outs, cpp_outs)

	s += output_diff_table_line("rastrigin_2d", py_outs, cpp_outs)
	s += output_diff_table_line("rastrigin_10d", py_outs, cpp_outs)
	s += output_diff_table_line("rastrigin_30d", py_outs, cpp_outs)

	s += output_diff_table_line("schwefel_2d", py_outs, cpp_outs)
	s += output_diff_table_line("schwefel_10d", py_outs, cpp_outs)
	s += output_diff_table_line("schwefel_30d", py_outs, cpp_outs)

	s += output_diff_table_line("rosenbrock", py_outs, cpp_outs)
	s += output_diff_table_line("schaffer_f6", py_outs, cpp_outs)

	return s


def output_table_line(function_name, cpp_outs):
	cpp_csv = list(filter(lambda x: function_name in x, cpp_outs))

	if len(cpp_csv) != 1:
		print("Cannot find csv file")
		return
	
	cpp_df = pd.read_csv("outputs/" + cpp_csv[0])

	s = ""

	s += "{}|{:.1f} \\textpm {:.1f} |{:.2e} \\textpm {:.2e}|{:.2e}\\textpm {:.2e}|{}\\textpm {}\n".format(
		function_name.replace("_", "\\_"),
		cpp_df["nit"].mean(), cpp_df["nit"].std(),
		cpp_df["x_diff"].mean(), cpp_df["x_diff"].std(),
		cpp_df["f_diff"].mean(), cpp_df["f_diff"].std(),
		duration_to_str(cpp_df["duration_microsec"].mean()), duration_to_str(cpp_df["duration_microsec"].std()),
		)
	
	return s

def output_table(method_name):
	files = os.listdir("outputs")
	outputs = list(filter(lambda x: ".csv" in x, files))
	cpp_outs = list(filter(lambda x: "cpp" in x and method_name in x, outputs))

	s = "Method|Iters|x\\_diff|f\\_diff|duration\n"
	s += "-|-|-|-|-\n"

	s += output_table_line("sphere_2d", cpp_outs)
	s += output_table_line("sphere_10d", cpp_outs)
	s += output_table_line("sphere_30d", cpp_outs)

	s += output_table_line("rastrigin_2d", cpp_outs)
	s += output_table_line("rastrigin_10d", cpp_outs)
	s += output_table_line("rastrigin_30d", cpp_outs)

	s += output_table_line("schwefel_2d", cpp_outs)
	s += output_table_line("schwefel_10d", cpp_outs)
	s += output_table_line("schwefel_30d", cpp_outs)

	s += output_table_line("rosenbrock", cpp_outs)
	s += output_table_line("schaffer_f6", cpp_outs)

	return s

if __name__ == "__main__":
	# box_plot_lowdim("L-BFGS-B", "lbfgs30", "L-BFGS-30")
	# box_plot_lowdim("Nelder-Mead", "nelder_mead", "Nelder-Mead")
	# box_plot_lowdim("CG", "cg_fr", "CG_FR")
	# box_plot_lowdim("CG", "cg_pr", "CG_PR")
	# box_plot_lowdim("CG", "cg_prp", "CG_PRP")

	s  = "### {}\n{}\n".format("LBFGS30", output_diff_table("L-BFGS-B", "lbfgs30") )
	s += "### {}\n{}\n".format("Nelder-Mead", output_diff_table("Nelder-Mead", "nelder_mead") )
	s += "### {}\n{}\n".format("CGFR", output_diff_table("CG", "cg_fr") )
	s += "### {}\n{}\n".format("CGPR", output_diff_table("CG", "cg_pr") )
	s += "### {}\n{}\n".format("CGPRP", output_diff_table("CG", "cg_prp") )
	open("figures/nm_tbl.md", "w").write(s)

	# s = output_table("genetic_bin_co") + "\n"
	# s += output_table("genetic_bin_gaus") + "\n"
	# s += output_table("genetic_rk_co") + "\n"
	# s += output_table("genetic_rk_gaus") + "\n"
	# open("figures/genetics.md", "w").write(s)

	# plot_initial_points(rosenbrock, 10, "rosenbrock")
	# plot_initial_points(schaffer_f6, 10, "schaffer_f6")
	# plot_initial_points(rastrigin, 5.12, "2_rastrigin")
	# plot_initial_points(sphere, 5, "2_sphere")
	# plot_initial_points(schwefel, 500, "2_schwefel")