#!/usr/bin/env python3

import animatplot as amp
import numpy as np
from matplotlib import pyplot as plt


def load_data():
    t = np.loadtxt("advection_t.dat")
    z = np.loadtxt("advection_z.dat")
    f = np.loadtxt("advection_f.dat")
    return t, z, f


def animate_f(t, z, f, save_as=None, fps=10):
    fig, ax = plt.subplots()

    line_block = amp.blocks.Line(z, f, ax=ax)

    timeline = amp.Timeline(t, fps=fps)

    animation = amp.Animation([line_block], timeline)

    ax.set_xlabel("z")
    ax.set_ylabel("f")

    ax.set_ylim([f.min(), f.max()])

    animation.controls(timeline_slider_args={"text": "t"})

    if save_as is not None:
        from matplotlib.animation import PillowWriter

        animation.save(save_as + ".gif", writer=PillowWriter(fps=fps))

    return animation


def plot_f0(t, f):
    plt.figure()
    plt.plot(t, f[:, 0])
    plt.xlabel("t")
    plt.ylabel("f(z=0)")


if __name__ == "__main__":
    t, z, f = load_data()
    animation = animate_f(t, z, f)
    # plot_f0(t, f)
    plt.show()
