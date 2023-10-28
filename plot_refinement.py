from plot_figure1 import Spatial, Temporal
import matplotlib.pyplot as plt
import os
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc("text.latex", preamble="\\usepackage{amsmath} \\usepackage{amsfonts}")

markers = ['d', 'o', 'x', 'v', '>', 's', '^']
colors = ['C0', 'C1', 'black', 'C2', 'C3', 'C4', 'C6']

def plot_spatial(ax):
    files = [
        os.path.join("out", "spatial-refinement-theta=0.5-kappa=-0.8-e=1.csv"),
        os.path.join("out", "spatial-refinement-theta=0.5-kappa=-0.4-e=1.csv"),
        os.path.join("out", "spatial-refinement-theta=0.5-kappa=0-e=1.csv"),
        os.path.join("out", "spatial-refinement-theta=0.5-kappa=1-e=1.csv"),
        os.path.join("out", "spatial-refinement-theta=0.5-kappa=8-e=1.csv"),
        os.path.join("out", "spatial-refinement-theta=0.5-kappa=64-e=1.csv"),
        os.path.join("out", "spatial-refinement-theta=0.5-kappa=256-e=1.csv")
    ]
    for f, m, c in zip(files, markers, colors):
        space = Spatial(f)
        point_kwargs = dict(marker=m, color=c, ls='None', mfc='None')
        line_kwargs = dict(ls='solid', color=c)
        if (abs(space.kappa - 8) < 1e-8) or (abs(space.kappa - 64) < 1e-8) or (abs(space.kappa - 256) < 1e-8):
            space.plot(ax, point_kwargs, line_kwargs, label_line=True)
        else:
            space.plot(ax, point_kwargs, line_kwargs, label_line=False)
        ax.annotate("$\\dfrac{M}{\\tau_{\\star}}=%i$" % ((space.m-1)/space.tau_star), xy=(0.6, 0.09), xycoords='axes fraction')


def plot_temporal(ax):
    files = [
        os.path.join("out", "temporal-refinement-theta=0.5-kappa=-0.8-e=1.csv"),
        os.path.join("out", "temporal-refinement-theta=0.5-kappa=-0.4-e=1.csv"),
        os.path.join("out", "temporal-refinement-theta=0.5-kappa=0-e=1.csv"),
        os.path.join("out", "temporal-refinement-theta=0.5-kappa=1-e=1.csv"),
        os.path.join("out", "temporal-refinement-theta=0.5-kappa=8-e=1.csv"),
        os.path.join("out", "temporal-refinement-theta=0.5-kappa=64-e=1.csv"),
        os.path.join("out", "temporal-refinement-theta=0.5-kappa=256-e=1.csv")
    ]
    for f, m, c in zip(files, markers, colors):
        I = Temporal(f)
        point_kwargs = dict(marker=m, color=c, ls='None', mfc='None', label='$\\kappa = %g$' % I.kappa)
        line_kwargs = dict(ls='solid', color=c)
        I.plot(ax, point_kwargs, line_kwargs, label_line=False)
        ax.annotate("$N=%i$" % ((I.n-1)), xy=(0.63, 0.05), xycoords='axes fraction')


if __name__ == '__main__':
    import os
    fig, ax = plt.subplots(ncols=2, figsize=(5.51181, 2.))
    plot_spatial(ax[0])
    plot_temporal(ax[1])
    ax[0].set_ylabel("$E_{\\chi,\\varphi}$", labelpad=0)
    ax[1].set_ylabel("$E^{\\tau,\\varphi}$", labelpad=0)
    leg = ax[0].legend(loc=(0.01, 0.62), frameon=False, edgecolor='None', facecolor='inherit', handlelength=0.7, handletextpad=0.2, fontsize="small")
    texts = list(leg.get_texts())
    texts[0].set_color(colors[4])
    texts[1].set_color(colors[5])
    texts[2].set_color(colors[6])
    ax[0].set_xlabel("$\\Delta \\chi = 1/N$", labelpad=0)
    ax[1].set_xlabel("$\\Delta \\tau = \\tau_{\\star}/M$", labelpad=0)
    leg = ax[1].legend(loc=(-1.5, 1.), ncols=7, handletextpad=0.1, columnspacing = 0.1, frameon=False, edgecolor='None', facecolor='inherit')
    for text, c in zip(leg.get_texts(), colors):
        text.set_color(c)
    for a in (ax[0], ax[1]):
        a.grid()
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    fig.subplots_adjust(bottom=0.2, right=0.99, wspace=0.28, left=0.11)
    fig.savefig('refinement.png', dpi=300, transparent=True)