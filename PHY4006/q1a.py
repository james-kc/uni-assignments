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
        df (DataFrame): Dataframe with columns 'photon_energy' and 'flux'.
    """

    _, ax = plt.subplots(figsize=(12, 6))

    ax.errorbar(
        x=df.photon_energy,
        y=df.flux,
        yerr=df.flux * .2,
        color='k',
        label='Cluster 3, 101 Mpc'
    )

    plt.legend()

    ax.set_yscale('log')
    ax.set_xlabel(
        'Photon Energy (keV)',
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

    plt.show()


def plot_luminosity(df):
    """Function for plotting X-Ray light curves from galexy clusters with flux
    measurement uncertainty of 20%.

    Args:
        df (DataFrame): Dataframe with columns 'photon_energy' and 'flux'.
    """

    _, ax = plt.subplots(figsize=(12, 6))

    ax.errorbar(
        x=df.photon_energy,
        y=df.luminosity,
        yerr=df.flux * 4 * const.pi * distance**2 * 0.2,
        color='k',
        label='Cluster 3, 101 Mpc'
    )

    plt.legend()

    ax.set_xlabel(
        'Photon Energy (keV)',
        fontsize=20
    )
    ax.set_ylabel(
        'Luminosity (J/s/Hz)',
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

    plt.show()


def plot_luminosity_fit(df):
    """Function for plotting X-Ray light curves from galexy clusters with flux
    measurement uncertainty of 20%.

    Args:
        df (DataFrame): Dataframe with columns 'photon_energy' and 'flux'.
    """

    _, ax = plt.subplots(figsize=(12, 6))

    ax.errorbar(
        x=df.frequency,
        y=np.log(df.luminosity),
        yerr=(df.flux**-1) * (0.2 * df.flux),
        color='k',
        label='Cluster 3, 101 Mpc'
    )

    ax.plot(
        df.frequency,
        df.fit,
        color='r',
        linestyle='--',
        label='Fit'
    )

    plt.legend()

    ax.set_xlabel(
        'Photon Frequency, $\\nu$ (Hz)',
        fontsize=20
    )
    ax.set_ylabel(
        'Log Luminosity (J/s/Hz)',
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

    plt.show()


def bremsstrahlung_luminosity(E, T):#, n_i, n_e):
    v = (E * 1e3 * const.e) / const.h
    Z = 26
    n_i = 1e10
    n_e = 1e10
    # T = 1e11

    constant = 6.8e-51
    Z_T_terms = Z**2 * T**(-1/2)
    density_terms = n_i * n_e
    exp_term = np.exp((-const.h * v) / (const.Boltzmann * T))
    
    return constant * Z_T_terms * density_terms * exp_term


def energy_to_wavelength(E):
    E_joules = E * 1e3 * const.e
    return (const.c * const.h) / E_joules


def power_law(E, T):

    wavelength = energy_to_wavelength(E)

    L_top = 2 * const.h * const.c**2

    L_bottom = wavelength**5

    R_bottom = np.exp(
        (const.h * const.c) /
        (wavelength * const.k * T)
    ) - 1

    curve = (L_top / L_bottom) * (1 / R_bottom)
    return curve / curve.max()


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

df['luminosity'] = df.flux * 4 * const.pi * distance**2
df['frequency'] = (df.photon_energy * 1e3 * const.e) / const.h

popt, pcov = curve_fit(
    straight_line,
    df.frequency,
    np.log(df.luminosity),
    sigma=(df.flux**-1) * (0.2 * df.flux)
)

print(popt)

grad = popt[0]
temperature = - const.h / (const.k * grad)

print(f"Temperature = {temperature:.2E} K")

df['fit'] = straight_line(df.frequency, *popt)

# plot_flux(df)
plot_luminosity_fit(df)
