from src.plotting_util import save_figure, Spatial, Temporal
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc("text.latex", preamble="\\usepackage{amsmath} \\usepackage{amsfonts}")


def plot_convergence(spatial_file, temporal_file1, temporal_file2, temporal_file3):
    fig, ax = plt.subplots(ncols=2, figsize=(6.5, 5.))

    space = Spatial(spatial_file)
    point_kwargs = dict(marker='o', mfc='None', color='C0', ls='None')
    point_kwargs['label'] = "$\\varepsilon=%g$" % space.e
    line_kwargs = dict(ls='solid', color='C0')
    space.plot(ax[0], point_kwargs, line_kwargs, label_line=True)
    ax[0].annotate("$M\\varepsilon=%g$" % ((space.m-1)*space.e), xy=(0.63, 0.05), xycoords='axes fraction')


    temp = Temporal(temporal_file1)
    point_kwargs = dict(marker='o', mfc='None', color='C0', ls='None')
    point_kwargs['label'] = "$\\varepsilon=%g$" % temp.e
    line_kwargs = dict(ls='solid', color='C0')
    ax[1].annotate("$N\\varepsilon=%g$" % ((temp.n-1)*temp.e), xy=(0.63, 0.05), xycoords='axes fraction')
    temp.plot(ax[1], point_kwargs, line_kwargs, label_line=True)

    temp = Temporal(temporal_file2)
    point_kwargs = dict(marker='x', mfc='None', color='C1', ls='None')
    point_kwargs['label'] = "$\\varepsilon=%g$" % temp.e
    line_kwargs = dict(ls='dashed', color='C1')
    ax[1].annotate("$N\\varepsilon=%g$" % ((temp.n-1)*temp.e), xy=(0.63, 0.05), xycoords='axes fraction')
    temp.plot(ax[1], point_kwargs, line_kwargs, label_line=False)

    temp = Temporal(temporal_file3)
    point_kwargs = dict(marker='d', mfc='None', color='C2', ls='None')
    point_kwargs['label'] = "$\\varepsilon=%g$" % temp.e
    line_kwargs = dict(ls='dotted', color='C2')
    ax[1].annotate("$N\\varepsilon=%g$" % ((temp.n-1)*temp.e), xy=(0.63, 0.05), xycoords='axes fraction')
    temp.plot(ax[1], point_kwargs, line_kwargs, label_line=False)

    ax[0].set_xlabel("$\\Delta \\chi/\\varepsilon=1/\\left(N\\varepsilon\\right)$", labelpad=0)
    ax[1].set_xlabel("$\\Delta \\tau / \\varepsilon=2/\\left(M\\varepsilon\\right)$", labelpad=0)
    ax[0].set_ylabel(r"$E\left(N, M; \varepsilon\right)= \max\limits_{\substack{i\in \mathcal{I}\\j\in\mathcal{J}}} \left\lvert a_i^j - W_i^j \right\rvert$")
    for a in (ax[0], ax[1]):
        a.grid()
        a.tick_params(which='both', direction='in')
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)

    return fig, ax


if __name__ == '__main__':
    fig, ax = plot_convergence(
        "out/spatial.csv", "out/temporal1.csv", 
        "out/temporal2.csv", "out/temporal3.csv", 
    )
    ax[0].legend(facecolor='None', frameon=False)
    ax[1].legend(loc=(0.03, 0.75), facecolor='None', frameon=False)
    fig.suptitle("Numerical experiments suggest\n second-order convergence to analytical solution")
    fig.subplots_adjust(bottom=0.1, right=0.97, top=0.9, wspace=0.2, left=0.14)
    for axes in ax:
        axes.set_xlim([0.01, 1])
    save_figure(fig, "figureS1.png")