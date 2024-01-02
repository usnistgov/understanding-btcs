from src.plotting_util import save_figure, TWConvergence, ISOTHERM
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc("text.latex", preamble="\\usepackage{amsmath} \\usepackage{amsfonts} \\usepackage{euscript}")


if __name__ == '__main__':
    fig = plt.figure(figsize=(5.512, 2.5), dpi=300)
    top = 0.82
    bot1 = 0.15
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
    
    axes[1].set_ylabel("$\\overline{\\EuScript{E}}_N$", rotation=0., labelpad=12)
    for i in range(tw.es.shape[0]):
        axes[0].annotate("$\\varepsilon = %3.2f$" % tw.es[i], 
        xy=(0.03, 0.4 + 0.1*i), xycoords="axes fraction", color=tw.colors[i])
    
    fig.suptitle("Numerical breakthrough concentrations approach boundary-layer theory\n as $\\varepsilon$ decreases, but only when $\\varepsilon$ is sufficiently large")
    save_figure(fig, "figure4.png")