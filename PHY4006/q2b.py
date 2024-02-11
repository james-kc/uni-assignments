import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from pprint import pprint



def flux_plotter(df, t_cutoff, nt_cutoff):
    _, ax = plt.subplots(figsize=(8,7))

    df.plot(
        ax=ax,
        x='energy',
        y='photon_flux'
    )

    ax.axvspan(
        df.energy.min(),
        t_cutoff,
        color='grey',
        alpha=0.5
    )

    ax.axvspan(
        nt_cutoff,
        df.energy.max(),
        color='grey',
        alpha=0.5
    )

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(
        'Energy (keV)',
        fontsize=20
    )
    ax.set_ylabel(
        'Flux (photons $s^{-1} cm^{-2} keV^{-1}$)',
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
        f'PHY4006/flux_curve_q2b.png',
        dpi=300,
        bbox_inches='tight'
    )

    return ax


def flux_fit_plotter(df, t_cutoff, nt_cutoff):
    ax = flux_plotter(df, t_cutoff, nt_cutoff)
    
    df.plot(
        ax=ax,
        x='energy',
        y='fit'
    )

    plt.savefig(
        f'PHY4006/flux_curve_fit_q2b.png',
        dpi=300,
        bbox_inches='tight'
    )

    return ax


def straight_line(x, m, c):
    return m * x + c


filename = 'PHY4006/13-May-2005_full_photon_data.csv'

df = pd.read_csv(filename)

# There is a gross 0 flux that gives -np.inf when np.log
# curve_fit() doesn't like this so just remove the data point
df = df[df.photon_flux != 0]

'''
Columns:

energy (keV)
photon_flux (ph/s/cm^2/keV)
'''

thermal_cutoff = 19
nonthermal_cutoff = 160

nonthermal_spectra_mask = (
    (df.energy >= thermal_cutoff) &
    (df.energy <= nonthermal_cutoff)
)

popt, pcov = curve_fit(
    straight_line,
    np.log(df.energy[nonthermal_spectra_mask]),
    np.log(df.photon_flux[nonthermal_spectra_mask])
)

print(popt)

df['fit'] = np.exp(
    straight_line(
        np.log(df.energy[nonthermal_spectra_mask]),
        *popt
    )
)

flux_fit_plotter(df, thermal_cutoff, nonthermal_cutoff)

plt.show()