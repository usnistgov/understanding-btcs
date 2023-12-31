from src.plotting_util import calculate_slope_error, save_figure
import matplotlib.pyplot as plt
import typing
from scipy.stats import linregress
import numpy as np
import pandas as pd

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc("text.latex", preamble="\\usepackage{amsmath} \\usepackage{amsfonts}")


class Convergence:
    def __init__(self, file: str):
        data = pd.read_csv(file)
        assert len(set(data['eps'])) == 1, "Not same epsilon"
        assert len(set(data['tau_star'])) == 1, "Not same tau_star"
        assert len(set(data['theta'])) == 1, "Not same theta"

        if 'kappa' in data.keys():
            param_key = 'kappa'
        elif 'isoparam0' in data.keys():
            param_key = 'isoparam0'
        else:
            print("key not found!")
        assert len(set(data[param_key])) == 1, "Not correct isotherm"
        self.kappa = data[param_key][0]
        self.e = data['eps'].to_numpy()[0]
        self.ns = data['n'].to_numpy()
        self.ms = data['m'].to_numpy()
        self.tau_star = data['tau_star'][0]
        self.du_inf = data['error_inf'].to_numpy()
    
    def plot(self, ax, symbol_kwargs: dict, line_kwargs: dict, label_line=False):
        if isinstance(self.ds, list):
            self.ds = np.array(self.ds)

        slope, intercept, r, p, stderr = linregress(np.log(self.ds), np.log(self.du_inf))
        slope_error = calculate_slope_error(np.log(self.ds), np.log(self.du_inf), slope, intercept)
        d_line = np.linspace(self.ds.min(), self.ds.max())
        label = None
        if label_line:
            label = "slope $\\approx %3.2f_{%3.2f}$" % (slope, slope_error)
            
        ax.loglog(d_line, d_line**slope * np.exp(intercept), label=label, **line_kwargs)

        # plot points second
        ax.loglog(self.ds, self.du_inf, **symbol_kwargs)


class Spatial(Convergence):
    def __init__(self, file: str):
        Convergence.__init__(self, file)
        self.ds = 1./(self.ns - 1)/self.e
        assert len(set(self.ms - 1)) == 1, "Not same ms {}".format(set(self.ms-1))
        self.m = self.ms[0]


class Temporal(Convergence):
    def __init__(self, file: str):
        Convergence.__init__(self, file)
        self.ds = self.tau_star/(self.ms - 1)/self.e
        assert len(set(self.ns - 1)) == 1, "Not same ns"
        self.n = self.ns[0]
    

def plot_convergence(spatial_file, temporal_file1, temporal_file2, temporal_file3):
    fig, ax = plt.subplots(ncols=2, figsize=(5.51181, 2.))

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
    # fig, ax = plot_convergence(
    #     "out/spatial.csv", "out/temporal1.csv", 
    #     "out/temporal2.csv", "out/temporal3.csv", 
    # )
    # for a in (ax[0], ax[1]):
    #     a.set_xlim([0.01, 1])

    # ax[0].legend(facecolor='None', frameon=False, fontsize='small')
    # ax[1].legend(loc=(0.03, 0.6), facecolor='None', frameon=False, fontsize='small')
    # fig.subplots_adjust(bottom=0.19, right=0.99, top=0.95, wspace=0.2, left=0.11)
    # fig.savefig("convergence-Nonlinear-0.5-0.5.png", transparent=True, dpi=300)

    fig, ax = plot_convergence(
        "out/spatial.csv", "out/temporal1.csv", 
        "out/temporal2.csv", "out/temporal3.csv", 
    )
    ax[0].legend(facecolor='None', frameon=False, fontsize='small')
    ax[1].legend(loc=(0.03, 0.6), facecolor='None', frameon=False, fontsize='small')
    fig.subplots_adjust(bottom=0.17, right=0.97, top=0.95, wspace=0.2, left=0.14)
    for axes in ax:
        axes.set_xlim([0.01, 1])
    save_figure(fig, "figureS1.png")