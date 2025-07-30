import numpy as np
import numpy.random as npr

DIR = "data/lowdim_inits/"
NUM_INIT = 30

def gen_schaffer_f6():
	matrix = npr.normal(size=(NUM_INIT, 2)) * 5
	np.savetxt(DIR + "schaffer_f6.csv", matrix, delimiter=',', fmt='%f')

def gen_rosenbrock():
	matrix = np.tile(np.array([-1.2, 1.0]), (NUM_INIT, 1)) + npr.normal(size=(NUM_INIT, 2)) * 0.5
	np.savetxt(DIR + "rosenbrock.csv", matrix, delimiter=',', fmt='%f')

def gen_sphere(dim: int):
	matrix = npr.normal(size=(NUM_INIT, dim)) * 1.0
	np.savetxt(DIR + f"{dim}_" + "sphere.csv", matrix, delimiter=',', fmt='%f')

def gen_schwefel(dim: int):
	matrix = npr.normal(size=(NUM_INIT, dim)) * 100.0
	np.savetxt(DIR + f"{dim}_" + "schwefel.csv", matrix, delimiter=',', fmt='%f')

def gen_rastrigin(dim: int):
	matrix = npr.normal(size=(NUM_INIT, dim)) * 2.0
	np.savetxt(DIR + f"{dim}_" + "rastrigin.csv", matrix, delimiter=',', fmt='%f')

if __name__ == "__main__":
	gen_rosenbrock()

	gen_sphere(2)
	gen_sphere(10)
	gen_sphere(30)

	gen_schwefel(2)
	gen_schwefel(10)
	gen_schwefel(30)

	gen_rastrigin(2)
	gen_rastrigin(10)
	gen_rastrigin(30)

	gen_schaffer_f6()