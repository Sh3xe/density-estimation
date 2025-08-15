import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as mcolors
import os

font_size = 14
cverr = pd.read_csv("outputs/st_hm_cverr.csv")["x"]
space = np.log10(pd.read_csv("outputs/st_hm_space.csv")["x"])
time = np.log10(pd.read_csv("outputs/st_hm_time.csv")["x"])

def plot_heatmap():
	data = cverr.values.reshape( (len(space), len(time)) )
	plt.figure(figsize=(10, 8))
	plt.imshow(
		data,
		alpha=0.7,
		cmap="viridis",
		interpolation="bilinear",
		aspect="auto",
		norm=LogNorm(vmin=cverr.min(), vmax=cverr.max()),
		extent=[min(space), max(space), min(time), max(time)]
	)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.colorbar(label="CV Error")

def make_title(data_path: str):
	names = {
		"nm": "Nelder-Mead",
		"gm": "Gaussian mut.",
		"rs": "Rank sel.",
		"bts": "Binary tournament",
		"cm": "Crossover mut.",
	}
	keywords = []
	for (key, val) in names.items():
		if key in data_path:
			keywords.append(val)
	return ", ".join(keywords)

def plot_alg_points(data_path: str):
	plot_heatmap()

	df = pd.read_csv(data_path)
	df = df[
		(df['x'] >= min(space)) & (df['x'] <= max(space)) &
		(df['y'] >= min(time)) & (df['y'] <= max(time))
	]
	sizes  = np.linspace(30, 60, len(df))
	alphas = np.linspace(0.8, 1.0, len(df))
	cmap = mcolors.LinearSegmentedColormap.from_list("black_to_red", ["black", "red"])
	norm = plt.Normalize(vmin=min(alphas), vmax=max(alphas))
	colors = cmap(norm(alphas))

	plt.scatter(df["x"], df["y"], color=colors, edgecolors="white", linewidths=1.5, s=sizes, alpha=alphas)
	plt.scatter([0], [-1.5], color="black", edgecolors="white", linewidths=1.5, s=60)
	plt.title(make_title(data_path))

def plot_points():
	plt.rcParams.update({"font.size": font_size})
	plt.rcParams["figure.dpi"] = 150
	for file in os.listdir("outputs"):
		if "cal_" in file:
			plot_alg_points(f"outputs/{file}")
			plt.tight_layout()
			plt.savefig(file + ".png", bbox_inches="tight")
			plt.close()

def plot_evol():
	plt.rcParams.update({"font.size": font_size})
	plt.rcParams["figure.dpi"] = 150
	for file in os.listdir("outputs"):
		if "cal_" in file:
			evolution = pd.read_csv(f"outputs/{file}")["value"]
			title = make_title(file)
			best_so_far = evolution.cummin()
			plt.plot(best_so_far, label=title)

	plt.xlabel("Iteration")
	plt.ylabel("Best value so far")
	plt.yscale("log")
	plt.legend(fontsize=10)
	plt.grid(True)
	plt.savefig("evol.png", bbox_inches="tight")


if __name__ == "__main__":	
	plot_evol()