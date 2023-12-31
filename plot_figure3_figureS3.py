from src.plotting_util import save_figure, AsymptoticConvergence
from src.isotherms import Langmuir
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc("text.latex", preamble="\\usepackage{amsmath} \\usepackage{amsfonts}")


def main():
    fig = plt.figure(figsize=(5.512, 2.5), dpi=300)
    top = 0.82
    bot1 = 0.15
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
    axes[1].set_ylabel("")
    for i in range(rarefaction.es.shape[0]):
        axes[0].annotate("$\\varepsilon = %3.2f$" % rarefaction.es[i], xy=(0.55, 0.05 + 0.1*i), xycoords="axes fraction", color=rarefaction.colors[i])
    
    axes[2].annotate("$\\kappa=1$", xy=(0.2, 0.8), xycoords="axes fraction", color='C1')
    axes[2].annotate("$\\kappa=-\\dfrac{1}{2}$", xy=(0.5, 0.4), xycoords="axes fraction", color='C0')

    fig.suptitle('With decreasing $\\varepsilon$, numerical breakthrough concentrations\n approach rarefaction wave, but not shock wave.')
    save_figure(fig, "figure3.png")


def si():
    fig, ax = plt.subplots(figsize=(5., 5.), dpi=300)

    for (marker, color, k) in [
        ("o", "C0", -0.1),
        ("x", "C1", -0.2),
        ("*", "C2", -0.4),
        ("d", "C3", -0.8)
    ]:
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
