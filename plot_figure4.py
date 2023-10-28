from src.isotherms import Langmuir
from src.plotting_util import save_figure
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})


def main():
    e = 1/6
    tau_star = 0.6
    m = 2**13 + 1
    n = 2**13 + 1
    x = np.array([
        i/(n - 1) for i in range(n)
    ])
    tau = np.array([
        j*tau_star/(m - 1) for j in range(m)
    ])
    iso3, iso2, iso1 = Langmuir(1), Langmuir(6), Langmuir(36)

    fig = plt.figure(figsize=(5.512, 4.))
    left1 = 0.105
    bottom = 0.11
    spacing = 0.03
    top = 0.98
    height = (top - bottom - 2*spacing) / 3

    isowidth = height
    ax_iso1 = fig.add_axes([left1, bottom, isowidth, height])
    ax_iso2 = fig.add_axes([left1, bottom + height + spacing, isowidth, height])
    ax_iso3 = fig.add_axes([left1, bottom + 2*height + 2*spacing, isowidth, height])
    c = np.linspace(0, 1)
    for ax, iso, i in ((ax_iso1, iso1, 0), (ax_iso2, iso2, 1), (ax_iso3, iso3, 2)):
        ax.plot(
            c, iso.f(c), '-', color="C%i" % i
        )
        ax.plot(c, c, '--', color='grey')
        ax.set_ylabel("$F=L\\left(c; %2.f\\right)$" % iso.k)
        if i == 0:
            ax.set_xlabel("$c$")
        else:
            plt.setp(ax.get_xticklabels(), visible=False)

    right = 0.85
    left2 = 0.1
    xstart = left1 + isowidth + left2
    moviewidth = right - xstart
    ax_movie1 = fig.add_axes([xstart, bottom, moviewidth, height])
    ax_movie2 = fig.add_axes([xstart, bottom + height + spacing, moviewidth, height])
    ax_movie3 = fig.add_axes([xstart, bottom + 2*height + 2*spacing, moviewidth, height])

    ax_legend = fig.add_axes([0.95, bottom, 0.04, top - bottom])
    ax_legend.tick_params(length=0)
    ax_legend.set_xticks([])
    ax_legend.set_xticklabels([])

    # cls1.simulate()
    # cls2.simulate()
    # cls3.simulate()

    # make colormap
    colors = plt.cm.cool(np.linspace(0., 1., m))

    colors_plotted = []
    T_plotted = []
    for ax, iso, i in ((ax_movie1, iso1, 0), (ax_movie2, iso2, 1), (ax_movie3, iso3, 2)):
        for j in range(0, m, m//10):
            W_j = np.loadtxt("movie-%2.1f/movie-%i.txt" % (iso.k, j), skiprows=2)
            ax.plot(x/e, W_j, color=colors[j])
            if i == 0:
                colors_plotted.append(colors[j])
                T_plotted.append(tau[j]/e)

        ax.tick_params(axis="y", which='both', direction='in')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel("$W\\left(x_i,\\tau_j\\right)$")
        ax.set_ylim([0., 1.])
        ax.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xlim([0., 1/e])
        ax.grid()
        if i == 0:
            ax.set_xlabel("$x_i/\\varepsilon$")
        else:
            plt.setp(ax.get_xticklabels(), visible=False)

    for k in range(len(T_plotted)): 
        ax_legend.plot([0, 1], [T_plotted[k], T_plotted[k]], '-', color=colors_plotted[k], clip_on=False)

    ax_legend.set_ylabel("$\\tau_j/\\varepsilon$")
    ax_legend.set_yticks(T_plotted)
    ax_legend.set_ylim([np.min(T_plotted), np.max(T_plotted)])
    ax_legend.spines['top'].set_visible(False)
    ax_legend.spines['bottom'].set_visible(False)
    ax_legend.spines['right'].set_visible(False)
    save_figure(fig, "figure4.png")


if __name__ == '__main__':
    main()