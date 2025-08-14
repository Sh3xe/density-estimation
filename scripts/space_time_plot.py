import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

cverr = pd.read_csv("outputs/st_hm_cverr.csv")["x"]
space = pd.read_csv("outputs/st_hm_space.csv")["x"]
time = pd.read_csv("outputs/st_hm_time.csv")["x"]

data = cverr.values.reshape( (len(space), len(time)) )

plt.figure(figsize=(10, 8))
plt.imshow(
	data,
	cmap="viridis",
	interpolation="bilinear",
	aspect="auto",
	norm=LogNorm(vmin=cverr.min(), vmax=cverr.max()),
	extent=[min(space), max(space), min(time), max(time)]
)
plt.xticks([])
plt.yticks([])
plt.colorbar(label="Intensity" )
plt.xlabel("lambda space")
plt.ylabel("lambda time")
plt.show()