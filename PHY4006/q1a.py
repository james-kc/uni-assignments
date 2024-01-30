import pandas as pd
from matplotlib import pyplot as plt


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 13})


def plot_flux(df):
    _, ax = plt.subplots(figsize=(12, 6))

    ax.errorbar(
        x=df.photon_energy,
        y=df.flux,
        yerr=df.flux * .2,
        color='k',
        label='Cluster 3, 101 Mpc'
    )

    plt.legend()

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
        f'PHY4006/flux_curve.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()


filename = 'PHY4006/cluster_3.csv'
distance = 101 # Mpc

df = pd.read_csv(filename)

# Columns:
# photon_energy - keV
# flux - J/m^2/s/Hz - Uncertainty = 20%

# flux is dominated by thermal bremsstrahlung
# Gaunt factor \aprox 1

plot_flux(df)