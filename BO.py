import numpy as np
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
import function
import os
import argparse


class BasicBayesianOptimizer:
    def __init__(self, bounds, noise=1e-6):
        self.bounds = np.array(bounds)
        self.X = []
        self.Y = []
        self.noise = noise
        self.kernel = Matern(nu=2.5)
        self.gp = GaussianProcessRegressor(kernel=self.kernel, alpha=noise**2, normalize_y=True)

    def update(self, x_new, y_new):
        """Add a new observation."""
        self.X.append(np.array(x_new))
        self.Y.append(y_new)
        self.gp.fit(np.array(self.X), np.array(self.Y))

    def suggest(self):
        """Suggest the next point to sample using Expected Improvement."""
        def ei(x):
            x = np.array(x).reshape(1, -1)
            mu, sigma = self.gp.predict(x, return_std=True)
            mu_sample_opt = np.min(self.Y)

            with np.errstate(divide='warn'):
                imp = mu_sample_opt - mu
                Z = imp / (sigma + 1e-9)
                ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
                return -ei[0]

        x0 = np.random.uniform(self.bounds[:, 0], self.bounds[:, 1])
        res = minimize(ei, x0=x0, bounds=self.bounds, method='L-BFGS-B')
        return res.x




if __name__ == "__main__":
    from scipy.stats import norm
    import matplotlib.pyplot as plt
    parser = argparse.ArgumentParser(description="A simple script with arguments.")

    parser.add_argument("--name", type=str, help="name of function", required=True)
    args = parser.parse_args()

    black_box_function = function.square_sin
    
    file = os.path.join("BO_data", args.name, "dimension.csv")
    file = open(file)

    dim = int(file.readline())

    bounds = []

    # read bounds
    for i in range(dim):
        line = file.readline()
        parts = line.split(",")
        bounds.append((float(parts[0]), float(parts[1])))

    file.close()

    # read already calculated values
    xs = []
    ys = []

    y_min = 1000000
    x_min = None
    best_trial = -1

    trial = 0
    
    file = os.path.join("BO_data", args.name, "history.csv")
    file = open(file)

    for line in file:
        parts = line.split(",")
        # + 1 for the files
        assert len(parts) == (dim + 1)
        x = [float(n) for n in parts[0:dim]]
        y = float(parts[-1])

        if (y < y_min):
            y_min = y
            x_min = x
            best_trial = trial

        xs.append(x)
        ys.append(y)
        trial += 1

    bo = BasicBayesianOptimizer(bounds=bounds)
    
    print("Completed run #{}".format(len(xs)))

    for i, y in enumerate(ys):
        bo.update(xs[i], y)

    next_x = [str(n) for n in bo.suggest()]

    print("Next: {}".format(next_x))

    output_file = os.path.join("BO_data", args.name, "next.csv")
    output_file = open(output_file, mode='w')
    output_file.write(",".join(next_x))
    output_file.close()

    print("Best: ({}, {}) at trial {}".format(x_min, y_min, best_trial))
