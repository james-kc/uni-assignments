import pandas as pd
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 13})

filename = 'PHY4003/stopping_power.csv'

df = pd.read_csv(filename)


def plot_single_energy(energy_table):
    fig, ax = plt.subplots()

    energy_table['energy_deposited_normalised'] = (
        energy_table['energy_deposited'] /
        energy_table['energy_deposited'].max()
    )

    energy_table.plot(
        ax=ax,
        x='r',
        y='energy_deposited_normalised',
        style='k'
    )

    # ax.set_yscale('log')
    ax.set_ylim(ymin=0.005)
    ax.set_xlim(xmin=0)
    ax.get_legend().remove()
    ax.set_xlabel(
        'Depth (cm)',
        fontsize=20
    )
    ax.set_ylabel(
        'Relative Dose',
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
        f'PHY4003/dose_depth_curve.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()


def plot_multi_energy(energy_table, tumor_start, tumor_end):
    fig, ax = plt.subplots()

    particle_energy_columns = (
        energy_table
        .filter(regex=r'^particle_energy')
        .columns.tolist()
    )

    energy_deposited_columns = (
        energy_table
        .filter(regex=r'^energy_deposited_\d+\.\d+_normalised$')
        .columns.tolist()
    )

    energy_deposited_columns.append('cumulative_energy_deposited_normalised')

    
    for col in energy_deposited_columns:

        if col == 'cumulative_energy_deposited_normalised':
            style_str = 'r'
        else:
            style_str = 'k'

        print(col)

        energy_table.plot(
            ax=ax,
            x='r',
            y=col,
            style=style_str
        )

    ax.axvspan(tumor_start, tumor_end, alpha=0.2)

    # ax.set_yscale('log')
    ax.set_ylim(ymin=0.005)
    ax.set_xlim(xmin=0)
    ax.get_legend().remove()
    ax.set_xlabel(
        'Depth (cm)',
        fontsize=20
    )
    ax.set_ylabel(
        'Relative Dose',
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
        f'PHY4003/dose_depth_curve_multi_energy.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()


def find_nearest_stopping_power(energy):

    closest_stopping_power = df.iloc[
        (df['kinetic_energy'] - energy).abs().idxmin()
    ]

    return closest_stopping_power['stopping_power']


def adipose_iteration(energy, dr):
    """Function to calculate the energy deposited by alpha particles in adipose
    tissue.

    Args:
        energy (float): Energy of incident alpha particles (MeV)
        dr (float): Spatial steps to take through the calculations (cm)

    Returns:
        DataFrame: Array containing distance, energy of alpha particle and the
        energy deposted during the time step.
    """
    
    energy_list = []
    r = dr # cm
    density = 0.92 # g/cm^3

    while True:

        energy_deposit_rate = (
            find_nearest_stopping_power(energy) * density
        ) # MeV/cm
        energy_deposited = energy_deposit_rate * dr # MeV

        if (energy - energy_deposited) > 0:
            energy = energy - energy_deposited # MeV

            energy_list.append([r, energy, energy_deposited])

            r += dr

        elif energy > 0:
            remaining_r = energy / energy_deposit_rate
            energy_deposited = energy_deposit_rate * remaining_r
            energy = energy - energy_deposited

            energy_list.append([r + remaining_r, energy, energy_deposited])

            r += dr

        else:
            energy_list += [
                [r, 0, 0],
                [r + dr, 0, 0],
                [r + 2*dr, 0, 0],
                [r + 3*dr, 0, 0],
                [r + 4*dr, 0, 0]
            ]

            break

    return pd.DataFrame(
        energy_list,
        columns=['r', 'particle_energy', 'energy_deposited']
    )


def multi_energy_adipose_interaction(energies, dr, intensities=None):

    if not intensities:
        intensities = [1 for i in range(len(energies))]

    energies.sort(reverse=True)

    # Initialising table with first energy in list
    energy = energies[0]
    out = adipose_iteration(energy, dr)
    out = out.rename(
        columns = {
            'particle_energy': f'particle_energy_{energy}',
            'energy_deposited': f'energy_deposited_{energy}'
        }
    )

    for energy in energies[1:]:

        current_energy_cols = adipose_iteration(energy, dr)[
            ['particle_energy', 'energy_deposited']
        ]

        out[f'particle_energy_{energy}'] = current_energy_cols['particle_energy']
        out[f'energy_deposited_{energy}'] = current_energy_cols['energy_deposited']

    out['cumulative_energy_deposited'] = out[
        [f'energy_deposited_{energy}' for energy in energies]
    ].sum(axis=1)

    for i, energy in enumerate(energies):

        out[f'energy_deposited_{energy}_normalised'] = (
            out[f'energy_deposited_{energy}'] * intensities[i] /
            out['cumulative_energy_deposited'].max()
        )

    out[f'cumulative_energy_deposited_normalised'] = (
            out['cumulative_energy_deposited'] /
            out['cumulative_energy_deposited'].max()
        )

    return out



# out = adipose_iteration(800, .002)
# out = adipose_iteration(800, 1)

# out = adipose_iteration(652, .1)
# out = adipose_iteration(671, .1)

# out = multi_energy_adipose_interaction([652, 671], .1)
out = multi_energy_adipose_interaction(
    np.linspace(652, 671, 5).tolist(),
    .1,
    intensities=[1, 100/85, 100/75, 100/55, 100/4]
)

out.to_csv('PHY4003/dose_depth.csv', index=False)

plot_multi_energy(out, 19.5, 20.5)
