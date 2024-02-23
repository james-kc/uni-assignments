import math
import pandas as pd
from matplotlib import pyplot as plt
import scipy.constants as const
import numpy as np
from scipy.optimize import curve_fit


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 13})


def plot_flux(df):
    """Function for plotting X-Ray light curves from galexy clusters with flux
    measurement uncertainty of 20%.

    Args:
        df (DataFrame): Dataframe with columns 'frequency' and 'flux'.

    Returns:
        ax (plt.Axes): Axes for matplotlib plot.
    """

    _, ax = plt.subplots(figsize=(12, 6))

    ax.errorbar(
        x=df.frequency,
        y=df.flux,
        yerr=df.flux * .2,
        color='k',
        label='Cluster 3, 101 Mpc'
    )

    plt.legend()

    ax.set_yscale('log')
    ax.set_xlabel(
        'Photon Frequency, $\\nu$ (Hz)',
        fontsize=20
    )
    ax.set_ylabel(
        'Flux (J/$m^2$/s/Hz)',
        fontsize=20
    )
    ax.tick_params(
        axis='x',
        labelsize=13,
        rotation=0,
        direction='in',
        top=True,
        labelbottom=True
    )
    ax.tick_params(
        axis='y',
        labelsize=13,
        which='both',
        direction='in',
        right=True,
        labelbottom=True
    )

    plt.savefig(
        f'PHY4006/flux_curve_q1a.png',
        dpi=300,
        bbox_inches='tight'
    )

    return ax


def plot_flux_fit(df):
    """Function for plotting X-Ray light curves from galexy clusters and a
    bremsstrahlung emission fit overplotted with flux measurement uncertainty of
    20%.

    Args:
        df (DataFrame): Dataframe with columns 'frequency', 'flux' and 'fit.

    Returns:
        ax (plt.Axes): Axes for matplotlib plot.
    """

    ax = plot_flux(df)

    ax.plot(
        df.frequency,
        df.fit,
        color='r',
        linestyle='--',
        label='Fit'
    )

    plt.legend()

    plt.savefig(
        f'PHY4006/flux_curve_fit_q1a.png',
        dpi=300,
        bbox_inches='tight'
    )

    return ax


def straight_line(x, m, c):
    return m * x + c


filename = 'PHY4006/cluster_3.csv'
distance = 101 * 1e6 * const.parsec # Mpc (cluster_3)

df = pd.read_csv(filename)

'''
# Columns:
# photon_energy - keV
# flux - J/m^2/s/Hz - Uncertainty = 20%

# flux is dominated by thermal bremsstrahlung
# Gaunt factor ~ 1
'''

df['frequency'] = (df.photon_energy * 1e3 * const.e) / const.h
df['luminosity'] = 4 * const.pi * distance**2 * df.flux

popt, pcov = curve_fit(
    straight_line,
    df.frequency,
    np.log(df.flux),
    sigma=np.abs(1/df.flux) * (0.2 * df.flux)
)

print(popt)

grad = popt[0]
temperature = - const.h / (const.k * grad)

print(f"Temperature = {temperature:.2E} K")

df['fit'] = np.exp(straight_line(df.frequency, *popt))

plot_flux_fit(df)
plt.show()

popt, pcov = curve_fit(
    straight_line,
    df.frequency,
    np.log(df.luminosity),
    sigma=(df.flux**-1) * (0.2 * df.flux)
)

print(popt)
