#!/usr/bin/env python3

import h5py
from matplotlib import pyplot as plt
import numpy as np
from statistics import mean, stdev


trendlines = True

# plot_func = plt.loglog
# plot_func = plt.semilogx
# plot_func = plt.semilogy
plot_func = plt.plot

n_max = 128


def plot_scan(result_file=None):
    if result_file is None:
        result_file = "scan.h5"

    if isinstance(result_file, str):
        result_file = [result_file]

    r_list = []
    for this_file in result_file:
        if this_file.split(".")[1] == "h5":
            r_list.append(h5py.File(this_file, "r"))
        elif this_file.split(".")[1] == "txt":
            model_group = {}
            with open(this_file, "r") as f:
                for line in f:
                    n, time = line.split()
                    n = n
                    time = float(time)
                    if n not in model_group:
                        model_group[n] = {"times": []}
                    model_group[n]["times"].append(time)
            # calculate average times
            for n_group in model_group.values():
                times = n_group["times"]
                n_group["avg_time"] = [mean(times)]
            r_list.append({"juliafft": model_group})
        else:
            raise ValueError(f"Unrecognised result file type {this_file}")

    plot_save_name = "_".join(x.split(".")[0] for x in result_file) + "_timing.pdf"

    fig = plt.figure()
    cycle = iter(plt.rcParams["axes.prop_cycle"])
    max_list = []
    for r in r_list:
        for model, model_group in r.items():
            n = [int(i) for i in model_group if int(i) <= n_max]
            n.sort()
            avg_time = np.array([model_group[str(i)]["avg_time"][0] for i in n])
            scatter = np.array([stdev(model_group[str(i)]["times"]) for i in n])
            n = np.array(n)
            color = next(cycle)["color"]
            plot_func(n, avg_time, label=model, color=color)
            # if plot_func is plt.plot:
            plt.fill_between(
                n, avg_time - scatter, avg_time + scatter, color=color, alpha=0.3
            )

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

    plt.xlabel("Nz")
    plt.ylabel("run time (s)")

    plt.legend()

    plt.savefig(plot_save_name)

    julia = None
    c = None
    for r in r_list:
        if "juliafft" in r:
            julia = r["juliafft"]
        if "chebyshevfft" in r:
            c = r["chebyshevfft"]
    if julia is not None and c is not None:
        plt.figure()
        n_julia = [int(i) for i in julia]
        n_c = [int(i) for i in c]
        n = [n for n in n_julia if n in n_c]
        ratio = [julia[str(i)]["avg_time"][0] / c[str(i)]["avg_time"][0] for i in n]
        plt.plot(n, ratio)
        plt.xlabel("Nz")
        plt.ylabel("time(Julia) / time(C++)")

    plt.show()


if __name__ == "__main__":
    from sys import argv

    if len(argv) == 0:
        result_file = None
    elif len(argv) == 1:
        result_file = argv[1]
    else:
        result_file = argv[1:]

    plot_scan(result_file)
