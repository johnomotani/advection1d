#!/usr/bin/env python3

import h5py
import multiprocessing
from pathlib import Path
from statistics import mean
import subprocess
import toml

from plot_performance_scan import plot_scan

n_values = list(range(8, 129)) + [256]
# Means we start longest jobs first, so should load-balance a bit better
n_values.reverse()
models = [
    "chebyshevmatrix",
    "chebyshevfft",
    "chebyshevfft_r2r",
    "chebyshevfft_r2c",
]

inputs = [(model, n) for n in n_values for model in models]

base_directory = Path().absolute()
advection1d_path = Path(__file__).parent.joinpath("..", "advection1d").absolute()

# Use half of cpus to ensure system is not full - attempt to make sure timing is not
# affected by time used for system processes, etc.
n_procs = multiprocessing.cpu_count() // 2

with open("input.toml") as f:
    base_options = toml.load(f)


def run_sim(run_dir):
    # run the simulation, return the run time
    this_run = subprocess.run(
        str(advection1d_path), check=True, stdout=subprocess.PIPE, cwd=run_dir
    )

    # get the timing
    stdout_lines = this_run.stdout.decode().split("\n")
    return float(stdout_lines[-2].split()[-1][:-1])


def time_case(inputs):
    model, n = inputs
    print(model, n, flush=True)

    # set up input
    p = Path(f"{model}_N{n}")
    try:
        p.mkdir()
    except FileExistsError:
        return model, n, None
    options = base_options.copy()
    options["spatial_type"] = model
    options["N"] = n
    with p.joinpath("input.toml").open("w") as f:
        toml.dump(options, f)

    # repeat runs to make timing more reliable
    return model, n, [run_sim(p) for i in range(5)]


results_file = h5py.File("scan.h5", "a")
for model in models:
    if model not in results_file:
        results_file.create_group(model)

with multiprocessing.Pool(n_procs) as pool:
    for model, n, times in pool.imap_unordered(time_case, inputs):
        if times is None:
            # case has already been run, so skip
            continue

        model_group = results_file[model]

        if str(n) in model_group:
            # directory was deleted, so re-running for some reason.
            del n_group[str(n)]
        n_group = model_group.create_group(str(n))
        n_group.create_dataset("times", data=times, dtype=float)
        n_group.create_dataset("avg_time", data=[mean(times)], dtype=float)
        results_file.flush()

results_file.close()

plot_scan()
