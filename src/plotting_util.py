import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as st
from scipy.stats import linregress
import os

BASE_DIR = os.path.abspath(os.path.dirname(__file__))


def calculate_slope_error(x, y, slope, intercept, conf=0.95):
    e = y - (intercept + slope * x)
    ssr = (e * e).sum()

    mean_squared_error = ssr / (len(e) - 2)

    xr = x - x.mean()
    m = st.t.ppf((1 + conf) / 2., len(e) - 1)

    return m * np.sqrt(mean_squared_error / (xr * xr).sum())


def plot_data(x, y, q, key, ax, indexes=None):
    # first, just plot raw data, then do calcs
    ax.loglog(x, q, 'x', label='_nolabel_', color='C1')
    ax.loglog(x, y, 'o', label='_nolabel_', markerfacecolor='None', color='C0')

    if indexes is None:
        indexes = np.array(list(range(len(x))))

    log_x, log_y, log_q = map(np.log, (x[indexes], y[indexes], q[indexes]))

    if key[-1] == 't':
        latex_var = '\\tau'
    else:
        latex_var = '\\eta'

    x_line = np.linspace(x[indexes].min(), x[indexes].max())
    slope, intercept, r, p, stderr = linregress(log_x, log_y)
    slope_error = calculate_slope_error(log_x, log_y, slope, intercept)
    ax.loglog(x_line, x_line ** slope * np.exp(intercept), '-',
              label='$E_{Y} = \\mathcal{O}\\left(\\left(\\Delta %s \\right)^{%3.2f \\pm %3.2f}\\right)$' % (
              latex_var, slope, slope_error), color='C0')

    x_line = np.linspace(x[indexes].min(), x[indexes].max())
    slope, intercept, r, p, stderr = linregress(log_x, log_q)
    slope_error = calculate_slope_error(log_x, log_q, slope, intercept)
    ax.loglog(x_line, x_line ** slope * np.exp(intercept), '--',
              label='$E_{Q} = \\mathcal{O}\\left(\\left(\\Delta %s \\right)^{%3.2f \\pm %3.2f}\\right)$' % (
                  latex_var, slope, slope_error),
              color='C1')


def plot_on_axis(df, ax, key):
    def get_xyq(df):
        return df[key].to_numpy(), df['Y-inf'].to_numpy(), df['Q-inf'].to_numpy()

    x, y, q = get_xyq(df)
    plot_data(x, y, q, key, ax, indexes=np.where(df['Fit'] == True))


def plot_convergence(f_space: str = None, f_time: str = None, f_out: str = 'convergence_analysis.png'):
    assert (f_space is not None) or (f_time is not None), "Need at least one to plot"
    import pandas as pd

    fig = plt.figure(figsize=(5.0, 4.5))
    anno_kwargs = dict(xy=(-0.15, 0.90), xycoords='axes fraction', ha='left', va='top')

    num_rows = 2
    if f_space is None or f_time is None:
        num_rows = 1
    if f_space is not None:
        df_space = pd.read_csv(f_space)
        ax_space = fig.add_subplot(num_rows, 1, 1)
        plot_on_axis(df_space, ax_space, 'dx')
        ax_space.set_ylabel(r"$\|E_i\|_\infty$", rotation=0., labelpad=10.)
        ax_space.set_xlabel(r"$\Delta x$")
        ax_space.grid()
        ax_space.legend(loc='lower right')
        ax_space.annotate("(a)", **anno_kwargs)
    if f_time is not None:
        df_time = pd.read_csv(f_time)
        ax_time = fig.add_subplot(num_rows, 1, num_rows)
        plot_on_axis(df_time, ax_time, 'dt')
        ax_time.set_ylabel(r"$\|E_i\|_\infty$", rotation=0., labelpad=10.)
        ax_time.set_xlabel(r"$\Delta \tau$")
        ax_time.grid()
        ax_time.legend(loc='lower right')
        ax_time.annotate("(b)", **anno_kwargs)

    if num_rows == 2:
        fig.subplots_adjust(bottom=0.13, right=0.99, top=0.99, hspace=0.4, left=0.15)
    else:
        fig.subplots_adjust(bottom=0.15, right=0.99, top=0.99, left=0.2)

    fig.savefig(f_out)


def save_figure(fig, name):
    """Saves figure out/ directory

    Parameters
    ----------
    fig : matplotlib.pyplot.figure
        instance of figure
    name : str
        name of figure to be saved
    """
    fig.savefig(os.path.join(BASE_DIR, "..", "out", name), transparent=True, dpi=300)
