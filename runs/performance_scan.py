#!/usr/bin/env python3

import h5py
from pathlib import Path
from statistics import mean
import subprocess
import toml

from plot_performance_scan import plot_scan

n_values = list(range(8, 129)) + [256]

base_directory = Path().absolute()
advection1d_path = Path(__file__).parent.joinpath("..", "advection1d").absolute()

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


def time_case(model, n):
    print(model, n, flush=True)

    # set up input
    p = Path(f"{model}_N{n}")
    p.mkdir()
    options = base_options.copy()
    options["spatial_type"] = model
    options["N"] = n
    with p.joinpath("input.toml").open("w") as f:
        toml.dump(options, f)

    # repeat runs to make timing more reliable
    return [run_sim(p) for i in range(5)]


results_file = h5py.File("scan.h5", "a")
for model in [
    "chebyshevmatrix",
    "chebyshevfft",
    "chebyshevfft_r2r",
    "chebyshevfft_r2c",
]:
    try:
        model_group = results_file[model]
    except KeyError:
        model_group = results_file.create_group(model)

    for n in n_values:
        try:
            times = time_case(model, n)
        except FileExistsError:
            # case has already been run, so skip
            continue

        if str(n) in model_group:
            # directory was deleted, so re-running for some reason.
            del n_group[str(n)]
        n_group = model_group.create_group(str(n))
        n_group.create_dataset("times", data=times, dtype=float)
        n_group.create_dataset("avg_time", data=[mean(times)], dtype=float)
        results_file.flush()

results_file.close()

plot_scan()
