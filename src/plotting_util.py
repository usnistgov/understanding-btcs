"""! @brief Utilities for making plots
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy import stats as st
from scipy.stats import linregress
import os
from src.isotherms import Langmuir


def get_u_TW_Langmuir(xi):
    """! get traveling wave solution

    @param xi (np.array) positions to evaluate concentrations in boundary layer

    @returns concentrations inside boundary layer
    """
    xi[xi < -10] = -10
    xi[xi > 10] = 10

    P = np.sqrt(2) - 1
    A = P*np.exp(2*xi)
    return 1 + A - np.sqrt(A*(2 + A))

ISOTHERM = Langmuir(1)
BASE_DIR = os.path.abspath(os.path.dirname(__file__))


def calculate_slope_error(x, y, slope, intercept, conf=0.95):
    """!
    @brief calculate the error in the slope estimated with linear regression

    @param x np array of x-values of plot
    @param y numpy array of y-values of plot
    @param slope float representing slope calculated
    @param intercept float representing intercept calculated
    @param conf float in (0, 1) representing confidence fraction. Defaults to 0.95
    """
    e = y - (intercept + slope * x)
    ssr = (e * e).sum()

    mean_squared_error = ssr / (len(e) - 2)

    xr = x - x.mean()
    m = st.t.ppf((1 + conf) / 2., len(e) - 1)

    return m * np.sqrt(mean_squared_error / (xr * xr).sum())


def save_figure(fig, name):
    """! @brief saves figure to *out/* directory

    @param fig (matplotlib.pyplot.figure) instance of figure to be saved
    @param name (str) name of figure to be saved
    """
    fig.savefig(os.path.join(BASE_DIR, "..", "out", name), transparent=True, dpi=300)


class Convergence:
    """! Class for performing convergence analysis
        
    @param kappa (float) value of Langmuir isotherm parameter \f$\kappa\f$
    @param e (float) value of \f$\varepsilon\f$
    @param ns (np.array) value of \f$N\f$ for each calculation
    @param ms (np.array) value of \f$M\f$ for each calculation
    @param tau_star (float) value of \f$\tau_{\star}\f$ 
    @param du_inf (np.array) infinity norm of error for each calculation

    """
    def __init__(self, file: str):
        """!Initialize the class

        @param file name of data file to analyze

        @return instance of the class
        """
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
        """!plot the figure

        @param ax matplotlib axis
        @param symbol_kwargs key-word arguments for symbols
        @param line_kwargs line keyword arguments
        @param label_line whether or not to label line

        @return None
        """
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


class Temporal(Convergence):
    """! @brief class for plotting temporal convergence 

    @param n (float or int?) number of \f$N\f$ (spatial intervals ), the same for each simulation
    @see Convergence

    """
    def __init__(self, file: str):
        """! @brief initialize the class

        @return instancs of class
        """
        Convergence.__init__(self, file)
        self.ds = self.tau_star/(self.ms - 1)/self.e
        assert len(set(self.ns - 1)) == 1, "Not same ns"
        self.n = self.ns[0]


class Spatial(Convergence):
    """! @brief class for plotting temporal convergence 

    @param m (float or int?) number of \f$M\f$ (time-step intervals), the same for each simulation
    @see Convergence

    """
    def __init__(self, file: str):
        """! @brief initialize the class

        @return instancs of class
        """
        Convergence.__init__(self, file)
        self.ds = 1./(self.ns - 1)/self.e
        assert len(set(self.ms - 1)) == 1, "Not same ms {}".format(set(self.ms-1))
        self.m = self.ms[0]


class AsymptoticConvergence:
    """! Class for plotting convergence in epsilon (not N or M)
    
    @param isotherm (Langmuir) type of isotherm
    @param klasses array of solutions to store
    @param es array of epsilons, hard-coded to 0.02, 0.04, 0.08, 0.16, 0.32, 0.64
    @param colors, list of colors to plot, hard-coded to plt.cm.viridis_r
    @param tau_star, total dimensionless time to simulate

    """
    def __init__(self, isotherm: Langmuir, file_prefix, tau_star):
        """! initialize the class

        @param isotherm isotherm to use
        @param file_prefix (str) prefix for looking at files
        @param tau_star (float) total dimensionless time to simulate

        @returns instance of class
        """
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

    def get_other_theory(self, ts):
        """! get btc associated with other theory

        @param ts (iterable) list of physical times to get btc for

        @returns array of concentrations (btc) according to list of times
        """
        return [
            self.isotherm.get_c_local_equilibrium(1., ts[j])
            for j in range(len(ts))
        ]

    def plot_btcs(self, ax):
        """! plot breakthrough curves of numerical simulation and regular perturbation theory
        
        @param ax matplotlib axis for plotting
        """

        for i in range(self.es.shape[0]-1, -1, -1):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            ax.plot(taus, btc, '-', label="$\\varepsilon = %3.2f$" % self.es[i], color=self.colors[i])

        taus = np.linspace(0., self.tau_star, 10000)
        u_asympt = self.get_other_theory(taus + 1)
        ax.plot(taus, u_asympt, '--', color='red', linewidth=0.8)
        ax.set_xlabel("$\\tau_j$")
        ax.set_ylabel("$W_{N}^j$ or $\\bar{C}_{N}^j$")

    def plot_error(self, ax, plot_kwargs, line_kwargs=None):
        """! plot errors as a function of \f$\varepsilon\f$

        @param ax matplotlib axis
        @param plot_kwargs kwargs for plotting
        @param line_kwargs kwargs for lines

        """
        btc_error = []
        for i in range(len(self.es)):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            u_asympt = self.get_other_theory(taus + 1.)
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


class TWConvergence(AsymptoticConvergence):
    """! Class for plotting traveling wave convergence in epsilon (not N or M)
    
    @param isotherm (Langmuir) type of isotherm
    @param klasses array of solutions to store
    @param es array of epsilons, hard-coded to 0.02, 0.04, 0.08, 0.16, 0.32, 0.64
    @param colors, list of colors to plot, hard-coded to plt.cm.viridis_r
    @param tau_star, total dimensionless time to simulate

    """
    def __init__(self, isotherm: Langmuir, file_prefix, tau_star):
        AsymptoticConvergence.__init__(self, isotherm, file_prefix, tau_star)

    def get_other_theory(self, t, eps):
        """! get btc associated with other theory

        @param t (np.array) list of physical times to get btc for
        @param eps (float) value of $\f\varepsilon\f$ for calculation

        @returns array of concentrations (btc) according to list of times
        """
        return get_u_TW_Langmuir((1. - ISOTHERM.shock_lambda()*t)/eps)

    def plot_btcs(self, ax):
        """! plot breakthrough curves 
        
        @param ax matplotlib axis
        @return none
        """
        taus = np.linspace(0., self.tau_star, 10000)
        for i in range(self.es.shape[0]-1, -1, -1):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            ax.plot(taus, btc, '-', label="$\\varepsilon = %3.2f$" % self.es[i],
                    color=self.colors[i])
            u_asympt = self.get_other_theory(taus+1., self.es[i])
            ax.plot(taus, u_asympt, '--', color=self.colors[i], linewidth=0.8)

        ax.set_xlabel("$\\tau_{j}$")
        ax.set_ylabel("$W_{N}^j$ or $\\bar{U}_{N}^j$")

    def plot_error(self, ax, plot_kwargs, line_kwargs=None):
        """! plot error
        @param ax matplotlib axis
        @param plot_kwargs kws for plotting
        @param line_kwargs kws for line, defaults to None
        @return none
        """
        btc_error = []
        for i in range(len(self.es)):
            btc = np.loadtxt(self.file_prefix + "%3.2f.dat" % self.es[i], skiprows=2)
            m = btc.shape[0]
            taus = np.array([j * self.tau_star/(m-1) for j in range(m)])
            u_asympt = self.get_other_theory(taus + 1, self.es[i])
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
