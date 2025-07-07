import json
import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio

def schwefel(x, y):
    s = 0.0
    for xi in [x, y]:
        if xi > 500 or xi < -500:
            s += 0.02 * xi * xi
        else:
            s += -xi * np.sin(np.sqrt(abs(xi)))
    return 418.9829 * 2 + s

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

def main():
    file_path = "population.json"
    output_gif_path = "outputs/population_animation.gif"
    data = load_population_data(file_path)
    images = create_population_plots(data)
    create_gif(images, output_gif_path)

if __name__ == "__main__":
    main()
