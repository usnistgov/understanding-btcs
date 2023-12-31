from plot_figure2_figureS1 import AsymptoticConvergence, calculate_slope_error
from src.plotting_util import save_figure
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

ISOTHERM = Langmuir(1)

def get_u_TW_Langmuir(xi):
    xi[xi < -10] = -10
    xi[xi > 10] = 10

    P = np.sqrt(2) - 1
    A = P*np.exp(2*xi)
    return 1 + A - np.sqrt(A*(2 + A))


class TWConvergence(AsymptoticConvergence):
    def __init__(self, isotherm: Langmuir, file_prefix, tau_star):
        AsymptoticConvergence.__init__(self, isotherm, file_prefix, tau_star)
    
    def get_other_theory(self, x, t, eps):
        return get_u_TW_Langmuir((x - ISOTHERM.shock_lambda()*t)/eps)

    def plot_btcs(self, ax):

        taus = np.linspace(0., self.tau_star, 10000)
        for i in range(self.es.shape[0]-1, -1, -1):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            ax.plot(taus, btc, '-', label="$\\varepsilon = %3.2f$" % self.es[i], 
                    color=self.colors[i])
            u_asympt = self.get_other_theory(1, taus+1, self.es[i])
            ax.plot(taus, u_asympt, '--', color=self.colors[i], linewidth=0.8)
        
        ax.set_xlabel("$\\tau_{j}$")
        ax.set_ylabel("$W_{N}^j$ or $\\bar{U}_{N}^j$")
    
    def plot_error(self, ax, plot_kwargs, line_kwargs=None):
        btc_error = []
        for i in range(len(self.es)):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            u_asympt = self.get_other_theory(1., taus + 1, self.es[i])
            error = np.max(np.abs(btc - u_asympt))
            btc_error.append(error)

        ax.semilogx(self.es, btc_error, **plot_kwargs)
        slope, intercept, r, p, stderr = linregress(np.log(self.es), np.log(btc_error))
        eps_line = np.linspace(self.es.min(), self.es.max())
        slope_error = calculate_slope_error(np.log(self.es), np.log(btc_error), slope, intercept)
        ax.grid()
        if line_kwargs is not None:
            label = "$a \\approx%3.2f_{%3.2f}$" % (slope, slope_error)
            if 'label' not in line_kwargs.keys():
                ax.annotate(label, xy=(0.05, 0.05), xycoords="axes fraction", color=line_kwargs['color'])
            else:
                line_kwargs['label'] += ', ' + label
            ax.loglog(eps_line, eps_line**slope * np.exp(intercept), **line_kwargs)

        ax.set_xlabel("$\\varepsilon$")
        ax.set_ylabel(r"$\max\limits_{j\in\mathcal{J}}\quad \left\lvert W_{N}^j - \bar{C}_{N}^j\right\rvert$")


if __name__ == '__main__':
    fig = plt.figure(figsize=(5.5, 2.), dpi=300)
    top = 0.96
    bot1 = 0.2
    width = 0.33
    axes = [
        fig.add_axes([0.1, bot1, 0.5 - 0.1, top-bot1]),
        fig.add_axes([0.45 + 0.18, bot1, 0.55 - 0.18-.01, top-bot1])
    ]

    tw = TWConvergence(ISOTHERM, "out/shock-", 2.)

    tw.plot_btcs(axes[0])
    tw.plot_error(
        axes[1], 
        dict(marker='d', color='black', ls='None', mfc='None', clip_on=False), 
    )
    for i in range(2):
        axes[i].tick_params(which='both', direction='in')
        axes[i].spines['top'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        axes[i].grid()
    
    
    # plt.setp(axes[1].get_yticklabels(), visible=False)
    # axes[0].legend(facecolor='inherit', frameon=False)
    axes[1].set_ylabel("$\\overline{\\mathcal{E}}_N$", rotation=0., labelpad=12)
    for i in range(tw.es.shape[0]):
        axes[0].annotate("$\\varepsilon = %3.2f$" % tw.es[i], 
        xy=(0.03, 0.4 + 0.1*i), xycoords="axes fraction", color=tw.colors[i])
    
    # fig.subplots_adjust(left=0.14, bottom=0.14, right=0.99, top=0.99)
    save_figure(fig, "figure4.png")