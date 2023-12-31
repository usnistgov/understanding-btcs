from src.plotting_util import calculate_slope_error, save_figure
import os
from src.isotherms import Langmuir
from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc("text.latex", preamble="\\usepackage{amsmath} \\usepackage{amsfonts}")


class AsymptoticConvergence:
    """
    t_hat : float
        total dimensionless time to simulate. Must be greater than 1
    """
    def __init__(self, isotherm: Langmuir, file_prefix, tau_star):
        self.isotherm = isotherm
        self.klasses = []  # storing solutions
        self.es = np.array([
            0.02, 0.04, 0.08, 0.16, 0.32,  0.64
        ])
        r = np.log10(self.es)
        r = (r - r.min())/(r.max() - r.min())
        self.colors = plt.cm.viridis_r(r)
        self.file_prefix = file_prefix
        self.tau_star = tau_star

    def get_other_theory(self, x, ts, eps):
        return [
            self.isotherm.get_c_local_equilibrium(1, ts[j])
            for j in range(len(ts))
        ]
    
    def plot_btcs(self, ax):

        for i in range(self.es.shape[0]-1, -1, -1):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            ax.plot(taus, btc, '-', label="$\\varepsilon = %3.2f$" % self.es[i], color=self.colors[i])
        
        taus = np.linspace(0., self.tau_star, 10000)
        u_asympt = self.get_other_theory(1, taus + 1, self.es.min())
        ax.plot(taus, u_asympt, '--', color='red', linewidth=0.8)
        ax.set_xlabel("$\\tau_j$")
        ax.set_ylabel("$W_{N}^j$ or $\\bar{C}_{N}^j$")
    
    def plot_error(self, ax, plot_kwargs, line_kwargs=None):
        btc_error = []
        for i in range(len(self.es)):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            u_asympt = self.get_other_theory(1., taus + 1, self.es[i])
            error = np.max(np.abs(btc - u_asympt))
            btc_error.append(error)

        ax.loglog(self.es, btc_error, **plot_kwargs)
        slope, intercept, r, p, stderr = linregress(np.log(self.es), np.log(btc_error))
        eps_line = np.linspace(self.es.min(), self.es.max())
        slope_error = calculate_slope_error(np.log(self.es), np.log(btc_error), slope, intercept)
        print(slope, slope_error)
        if line_kwargs is not None:
            label = "$a \\approx%3.2f_{%3.2f}$" % (slope, slope_error)
            if 'label' not in line_kwargs.keys():
                ax.annotate(label, xy=(0.05, 0.05), xycoords="axes fraction", color=line_kwargs['color'])
            else:
                line_kwargs['label'] += ', ' + label
            ax.loglog(eps_line, eps_line**slope * np.exp(intercept), **line_kwargs)

        ax.set_xlabel("$\\varepsilon$")
        ax.set_ylabel(r"$\max\limits_{j\in\mathcal{J}}\quad \left\lvert W_{N}^j - \bar{C}_{N}^j\right\rvert$")


def main():

    fig = plt.figure(figsize=(5.512, 2.), dpi=300)
    top = 0.97
    bot1 = 0.2
    left = 0.1
    space=0.02
    width = 0.33
    axes = [
        fig.add_axes([left, bot1, (2*width - left - space)/2, top-bot1]),
        fig.add_axes([left +(2*width - left)/2, bot1, (2*width - left - space)/2, top-bot1]),
        fig.add_axes([2*width + left, bot1, width - left, top-bot1])
    ]

    rarefaction = AsymptoticConvergence(Langmuir(-1/2), "out/rarefaction-", 4.)
    shock = AsymptoticConvergence(Langmuir(1), "out/shock-", 2.)

    rarefaction.plot_btcs(axes[0])
    shock.plot_btcs(axes[1])
    rarefaction.plot_error(
        axes[2], 
        dict(marker='o', color='C0', ls='None', mfc='None', clip_on=False), 
        dict(ls='dashed', color='C0')
    )
    shock.plot_error(
        axes[2], 
        dict(marker='x', color='C1', ls='None', clip_on=False)
    )
    for i in range(3):
        axes[i].tick_params(which='both', direction='in')
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        axes[i].grid()
    
    ticks = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
    axes[2].set_ylim([ticks[0], ticks[-1]])
    axes[2].set_yticks(ticks)
    axes[2].set_yticklabels(list(map(str, ticks)))
    for i in range(2):
        axes[i].set_ylim([-0.05, 1.05])
    
    plt.setp(axes[1].get_yticklabels(), visible=False)
    # axes[0].legend(facecolor='inherit', frameon=False)
    axes[1].set_ylabel("")
    for i in range(rarefaction.es.shape[0]):
        axes[0].annotate("$\\varepsilon = %3.2f$" % rarefaction.es[i], xy=(0.55, 0.05 + 0.1*i), xycoords="axes fraction", color=rarefaction.colors[i])
    
    axes[2].annotate("$\\kappa=1$", xy=(0.2, 0.8), xycoords="axes fraction", color='C1')
    axes[2].annotate("$\\kappa=-\\dfrac{1}{2}$", xy=(0.5, 0.4), xycoords="axes fraction", color='C0')

    # fig.subplots_adjust(left=0.08, bottom=0.14, right=0.99, top=0.99)
    save_figure(fig, "figure3.png")


def si():
    fig, ax = plt.subplots(figsize=(5., 5.), dpi=300)
    kwargs = dict(DX=0.01, DT=0.004, theta=0.5)

    for (marker, color, k) in [
        ("o", "C0", -0.1),
        ("x", "C1", -0.2),
        ("*", "C2", -0.4),
        ("d", "C3", -0.8)
    ]:
        isotherm = Langmuir(k)
        rarefaction = AsymptoticConvergence(Langmuir(k), "out/rarefaction-%2.1f-" % k, 4.)
        rarefaction.plot_error(
            ax, 
            dict(marker=marker, color=color, ls='None', mfc='None', clip_on=False), 
            dict(ls='solid', color=color, label="$\\kappa=%2.1f$" % k)
        )

    ax.legend(facecolor='inherit', frameon=False)
    ax.set_xlim([0.01, 1.])
    ax.set_ylim([0.1, 1])
    ax.grid()
    fig.subplots_adjust(left=0.2, top=0.97, right=0.97)
    save_figure(fig, "figureS3.png")




if __name__ == '__main__':
    main()
    si()
