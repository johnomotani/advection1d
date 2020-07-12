#!/usr/bin/env python3

import animatplot as amp
import numpy as np
from matplotlib import pyplot as plt


def load_data():
    t = np.loadtxt("advection_t.dat")
    f = np.loadtxt("advection_f.dat")
    return t, f


def animate_f(t, f, save_as=None, fps=10):
    fig, ax = plt.subplots()
    line_block = amp.blocks.Line(f, ax=ax)
    timeline = amp.Timeline(t, fps=fps)
    animation = amp.Animation([line_block], timeline)
    ax.set_xlabel("z-index")
    ax.set_ylabel("f")
    animation.controls(timeline_slider_args={"text": "t"})
    if save_as is not None:
        from matplotlib.animation import PillowWriter

        animation.save(save_as + ".gif", writer=PillowWriter(fps=fps))

    return animation


if __name__ == "__main__":
    t, f = load_data()
    animation = animate_f(t, f)
    plt.show()
