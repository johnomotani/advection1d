#!/usr/bin/env python3

import animatplot as amp
import numpy as np
from matplotlib import pyplot as plt
from sys import argv


def load_data_Michael(inputfilename):
    def commastrip(s):
        if s.endswith(b","):
            return s[:-1]
        else:
            return s

    t, z, v, f = np.loadtxt(
        inputfilename,
        usecols=(1, 3, 5, 7),
        converters={1: commastrip, 3: commastrip, 5: commastrip},
        unpack=True,
    )

    unique_ts = np.unique(t)
    unique_zs = np.unique(z)
    unique_vs = np.unique(v)

    nt = len(unique_ts)
    nz = len(unique_zs)
    nv = len(unique_vs)
    f = f.reshape(nt, nv, nz)

    # Only keep largest v for now in comparison to John's code...
    f = f[:, -1, :]

    return unique_ts, unique_zs, f


def load_data_John(directoryname):
    t = np.loadtxt(directoryname + "/advection_t.dat")
    z = np.loadtxt(directoryname + "/advection_z.dat")
    f = np.loadtxt(directoryname + "/advection_f.dat")
    return t, z, f


def make_line(z, f, label):
    line_block = amp.blocks.Line(z, f, label=label)

    return line_block


def animate_lines(line_blocks, t, ymin, ymax, save_as=None, fps=10):
    timeline = amp.Timeline(t, fps=fps)
    animation = amp.Animation(line_blocks, timeline)

    plt.xlabel("z")
    plt.ylabel("f")
    plt.ylim([ymin, ymax])

    legend = plt.legend(loc="upper left")
    legend.set_draggable(True)

    animation.controls(timeline_slider_args={"text": "t"})

    if save_as is not None:
        from matplotlib.animation import PillowWriter

        animation.save(save_as + ".gif", writer=PillowWriter(fps=fps))

    return animation


if __name__ == "__main__":
    # Load and plot data from Michael's 'moment_kinetics'
    #####################################################
    inputfilename_M = argv[1]
    t_M, z_M, f_M = load_data_Michael(inputfilename_M)
    line_M = make_line(z_M, f_M, "Michael")

    # Load and plot data from John's 'advection_1d'
    ###############################################
    inputdirectoryname_J = argv[2]
    t_J, z_J, f_J = load_data_John(inputdirectoryname_J)

    # offset z to make it consistent with Michael's
    z_J = z_J - z_J[0] + z_M[0]

    line_J = make_line(z_J, f_J, "John")

    # check times are consistent for animation
    if not np.allclose(t_M, t_J, rtol=1.0e-14):
        raise ValueError(f"times are not consistent\n{t_M}\n{t_J}")

    # get limits for plot
    fmin = max(f_M.min(), f_J.min())
    fmax = max(f_M.max(), f_J.max())

    anim = animate_lines([line_M, line_J], t_M, fmin, fmax, save_as="compare")

    plt.show()
