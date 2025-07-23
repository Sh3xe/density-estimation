import numpy as np
import matplotlib.pyplot as plt

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

def plot():
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

if __name__ == "__main__":
	plot()