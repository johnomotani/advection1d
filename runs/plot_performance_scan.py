#!/usr/bin/env python3

import h5py
from matplotlib import pyplot as plt
import numpy as np


trendlines = True

# plot_func = plt.loglog
# plot_func = plt.semilogx
# plot_func = plt.semilogy
plot_func = plt.plot

n_max = 128


def plot_scan(result_file=None):
    if result_file is None:
        result_file = "scan.h5"
    r = h5py.File(result_file, "r")

    plot_save_name = result_file.split(".")[0] + "_timing.pdf"

    fig = plt.figure()
    cycle = iter(plt.rcParams["axes.prop_cycle"])
    max_list = []
    for model, model_group in r.items():
        n = [int(i) for i in model_group if int(i) <= n_max]
        n.sort()
        avg_time = [model_group[str(i)]["avg_time"][0] for i in n]
        n = np.array(n)
        color = next(cycle)["color"]
        plot_func(n, avg_time, label=model, color=color)

        if model != "chebyshevmatrix":
            # exclude matrix solve from setting plot limits
            max_list.append(max(avg_time))

        if trendlines:
            # plot power law trend line
            if model == "chebyshevmatrix":
                fit_func = lambda x: x ** 2
            elif "chebyshevfft" in model:
                # fit_func = lambda x: x*np.log(x)
                fit_func = lambda x: x
            else:
                raise ValueError(f"Unrecognised model {model}")

            # find point that performs best compared to expected trend
            fit_values = fit_func(n.astype(float))
            best_ind = np.argmin(avg_time / fit_values)
            fit_values *= avg_time[best_ind] / fit_values[best_ind]

            plot_func(n, fit_values, "--", color=color)

    ymin, ymax = plt.ylim()
    plt.ylim([ymin, max(max_list)])

    plt.legend()

    plt.savefig(plot_save_name)

    plt.show()


if __name__ == "__main__":
    from sys import argv

    if len(argv) > 1:
        result_file = argv[1]
    else:
        result_file = None

    plot_scan(result_file)
