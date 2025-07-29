import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filename = "outputs/lbfgs30_snp500_log_density.csv"
log_dens = pd.read_csv(filename)
plt.plot(log_dens)
# plt.hist(pd.read_csv("data/snp500/data_space.csv")["x"], density=True, bins=30)
plt.show()