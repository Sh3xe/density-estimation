import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import pandas as pd
import json
import imageio
from matplotlib.patches import Patch

FDAPDE_COL = "#EF8A87"
SCIPY_COL = "#85A2EB"

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

def plot_function_graph(function, bound, png_name):
	x = np.arange(-bound, bound, bound/30.0)
	y = np.arange(-bound, bound, bound/30.0)
	X, Y = np.meshgrid(x, y)
	Z = np.vectorize(function)(X, Y)

	fig = plt.figure(figsize=(6, 6))
	ax = fig.add_subplot(111, projection='3d')
	surf = ax.plot_surface(X, Y, Z, cmap='viridis', rstride=1, cstride=1, edgecolor='none', shade=True)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.view_init(elev=45, azim=45)
	ax.set_proj_type('ortho')
	# plt.show()
	plt.savefig(f"figures/plot_{png_name}.png")
	plt.close()

def box_plot_data(axis, use_log, data, tick_labels, img_name):
	axis.set_ylabel("||f-f*||")
	medianprops = dict(linestyle='-', linewidth=1.0, color='black')
	if use_log:
		axis.set_yscale("log")

	for i in range(len(tick_labels)):
		tick_labels[i] = tick_labels[i].replace("_", " ")
	bplot = axis.boxplot( data, patch_artist=True, notch = False, tick_labels=tick_labels, widths=0.6, medianprops=medianprops )

	for i in range(len(data)//2):
		bplot["boxes"][2*i].set_facecolor(FDAPDE_COL) # Col for fdaPDE
		bplot["boxes"][2*i+1].set_facecolor(SCIPY_COL) # Col for SciPy

	legend_elements = [
		Patch(facecolor=FDAPDE_COL, label='fdaPDE'),
    Patch(facecolor=SCIPY_COL, label='SciPy')
	]

	axis.legend(handles=legend_elements)
	plt.tight_layout()
	plt.xticks(rotation=70)
	plt.savefig(f"figures/{img_name}.png", bbox_inches="tight")
	plt.close()

def box_plot_lowdim(py_name, cpp_name, fig_title):
	files = os.listdir("outputs")
	outputs = list(filter(lambda x: ".csv" in x, files))
	py_outs = list(filter(lambda x: "py" in x and py_name + "_" in x, outputs))
	cpp_outs = list(filter(lambda x: "cpp" in x and cpp_name  + "_" in x, outputs))

	# Carefull when changing the order !!!
	test_cases = ["sphere", "rastrigin", "schwefel"]
	dimensions = [2, 10, 30]
	other_test_cases = ["rosenbrock", "schaffer_f6"]
	use_logs = [True, True, True]

	for (test_case, use_log) in zip(test_cases, use_logs):
		data = []; tick_labels= []
		for i, dim in enumerate(dimensions):
			method_title = f"{test_case}_{dim}d"
			py_csv = list(filter(lambda x: method_title in x, py_outs))
			cpp_csv = list(filter(lambda x: method_title in x, cpp_outs))
			assert(len(py_csv) == 1 and len(cpp_csv) == 1)
			data.append( pd.read_csv("outputs/" + cpp_csv[0])["f_diff"].to_numpy() )
			tick_labels.append(f"{dim}D")
			data.append( pd.read_csv("outputs/" + py_csv[0])["f_diff"].to_numpy() )
			tick_labels.append(f"{dim}D")
			fig, axis = plt.subplots()
		box_plot_data(axis, use_log, data, tick_labels, f"lowdim_comp_{test_case}_{cpp_name}")
		plt.close()	

	fig, ax = plt.subplots(figsize=(6, 6))
	data = []
	for test_case in other_test_cases:
		py_csv = list(filter(lambda x: test_case in x, py_outs))
		cpp_csv = list(filter(lambda x: test_case in x, cpp_outs))
		assert(len(py_csv) == 1 and len(cpp_csv) == 1)
		data.append( pd.read_csv("outputs/" + cpp_csv[0])["f_diff"].to_numpy() )
		data.append( pd.read_csv("outputs/" + py_csv[0])["f_diff"].to_numpy() )
		
	box_plot_data(ax, True, data, ["Rs", "Rs", "Sc", "Sc"], "lowdim_comp_{}_{}".format("_".join(other_test_cases), cpp_name))
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
	plt.scatter(population["x"], population["y"], color="black", edgecolors="white", linewidths=1.5, s=30, alpha=0.7)
	plt.title(f"Population at Time Frame {time_frame}")
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.xlim(-MAX_VAL, MAX_VAL)
	plt.ylim(-MAX_VAL, MAX_VAL)
	plt.grid(True)

def create_population_plots(iter_list):
	images = []
	for time_frame in iter_list:
		plt.figure()
		population = pd.read_csv(f"genetic_csv/{time_frame}_pop.csv")
		plot_population(population, time_frame)
		plt.savefig(f"population_{time_frame}.png")
		images.append(imageio.imread(f"population_{time_frame}.png"))
	plt.close()
	return images

def create_gif(images, output_path):
	imageio.mimsave(output_path, images, fps=10)

def animation_genetic():
	output_gif_path = "outputs/population_animation.gif"
	images = create_population_plots([2, 10, 20])
	create_gif(images, output_gif_path)

def plot_initial_points(function, bound, csv_name):
	# Find initial points
	df = pd.read_csv(f"data/lowdim_inits/{csv_name}.csv", names=["x", "y"])

	# Plots function
	x = np.linspace(-bound, bound, 200)
	y = np.linspace(-bound, bound, 200)
	Z = np.array([[function(xi, yi) for xi in x] for yi in y])

	plt.imshow(Z, extent=[-bound, bound, -bound, bound], origin="lower", cmap="viridis", alpha=0.7)
	plt.scatter(x=df["x"], y=df["y"], color="white")
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.xlim(-bound, bound)
	plt.ylim(-bound, bound)
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

	s += "{} |fdaPDE|{:.1f} \\textpm {:.1f} |{:.2e}\\textpm {:.2e}|{}\\textpm {}\n".format(
		function_name.replace("_", "\\_"),
		cpp_df["nit"].mean(), cpp_df["nit"].std(),
		cpp_df["f_diff"].mean(), cpp_df["f_diff"].std(),
		duration_to_str(cpp_df["duration_microsec"].mean()), duration_to_str(cpp_df["duration_microsec"].std()),
		)
	
	s += "{} | scipy|{:.1f} \\textpm {:.1f} |{:.2e}\\textpm {:.2e}|{}\\textpm {}\n".format(
		function_name.replace("_", "\\_"),
		py_df["nit"].mean(), py_df["nit"].std(),
		py_df["f_diff"].mean(), py_df["f_diff"].std(),
		duration_to_str(py_df["duration_microsec"].mean()), duration_to_str(py_df["duration_microsec"].std()),
		)
	return s

def output_diff_table(py_name, cpp_name):
	files = os.listdir("outputs")
	outputs = list(filter(lambda x: ".csv" in x, files))
	py_outs = list(filter(lambda x: "py" in x and py_name + "_" in x, outputs))
	cpp_outs = list(filter(lambda x: "cpp" in x and cpp_name + "_" in x, outputs))

	s = "Method|Type|Iters|f\\_diff|duration\n"
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

def box_plot_comp_all_cpp(function_name, methods, plot_title, use_log=False):
	files = os.listdir("outputs")
	cpp_csvs = filter(lambda x: "cpp" in x and function_name in x, files)
	data = []
	for csv_file_name in cpp_csvs:
		for m in methods:
			if m + "_" in csv_file_name:
				data.append(pd.read_csv(f"outputs/{csv_file_name}")["x_diff"].to_numpy())
				break
	

	plt.figure(figsize=(14, 6))
	plt.ylabel("x_diff")
	if use_log:
		plt.yscale("log")

	bplot = plt.boxplot(
		data,
		patch_artist=True,
		tick_labels=methods
	)

	bplot["boxes"][1].set_facecolor("#85A2EB")
	plt.title(plot_title)
	plt.tight_layout()
	plt.savefig(f"figures/{function_name}_all_comparison.png", dpi=300)
	plt.close()

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
	plt.rcParams.update({'font.size': 14})
	plt.rcParams['figure.dpi'] = 300
	# methods = ["lbfgs30", "genetic_bin_co", "genetic_bin_gaus", "genetic_rk_gaus", "genetic_rk_co", "nelder_mead"]
	# box_plot_comp_all_cpp("rastrigin_10d", methods, "Rastrigin 10D")
	# box_plot_comp_all_cpp("schaffer_f6", methods, "Schaffer F6", True)
	# box_plot_comp_all_cpp("schwefel_10d", methods, "Schewefel 10D")

	# box_plot_lowdim("L-BFGS-B", "lbfgs30", "L-BFGS-30")
	# box_plot_lowdim("Nelder-Mead", "nelder_mead", "Nelder-Mead")
	# box_plot_lowdim("CG", "cg_pr_restart", "CG_PR")
	# box_plot_lowdim("CG", "cg_prp_restart", "CG_PRP")

	# s  = "### {}\n{}\n".format("LBFGS30", output_diff_table("L-BFGS-B", "lbfgs30") )
	# s += "### {}\n{}\n".format("Nelder-Mead", output_diff_table("Nelder-Mead", "nelder_mead") )
	# s += "### {}\n{}\n".format("CGPRP", output_diff_table("CG", "cg_pr_restart") )
	# open("figures/nm_tbl.md", "w").write(s)

	# s = output_table("genetic_bin_co") + "\n"
	# s += output_table("genetic_bin_gaus") + "\n"
	# s += output_table("genetic_rk_co") + "\n"
	# s += output_table("genetic_rk_gaus") + "\n"
	# open("figures/genetics.md", "w").write(s)

	# plot_initial_points(rosenbrock, 3, "rosenbrock")
	# plot_initial_points(schaffer_f6, 10, "schaffer_f6")
	# plot_initial_points(rastrigin, 5.12, "2_rastrigin")
	# plot_initial_points(sphere, 5, "2_sphere")
	# plot_initial_points(schwefel, 500, "2_schwefel")

	# plot_function_graph(rosenbrock, 3, "rosenbrock")
	# plot_function_graph(schaffer_f6, 10, "schaffer_f6")
	# plot_function_graph(rastrigin, 5.12, "2_rastrigin")
	# plot_function_graph(sphere, 5, "2_sphere")
	# plot_function_graph(schwefel, 500, "2_schwefel")

	animation_genetic()