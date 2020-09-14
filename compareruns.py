#!/usr/bin/env python3

import animatplot as amp
import numpy as np
from matplotlib import pyplot as plt
from sys import argv


def load_data_Michael(inputfilename):
    def commastrip(s):
        if s.endswith(b','):
            return s[:-1]
        else:
            return s

    t,z,v,f = np.loadtxt(
        inputfilename,
        usecols=(1,3,5,7),
        converters={1:commastrip, 3:commastrip, 5:commastrip},
        unpack=True
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


def load_data(directoryname):
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
    times = []
    fs = []
    lines = []
    for d in argv[1:]:
        t, z, f = load_data(d)

        times.append(t)
        fs.append(f)
        lines.append(make_line(z, f, d))

    # check times are consistent for animation
    for t in times[1:]:
        if not np.allclose(times[0], t, rtol=1.e-14):
            raise ValueError(f"times are not consistent\n{times[0]}\n{t}")

    # get limits for plot
    fmin = max(f.min() for f in fs)
    fmax = max(f.max() for f in fs)

    anim = animate_lines(lines, times[0], fmin, fmax)

    plt.show()
