import pandas as pd
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt

filename = 'PHY4003/stopping_power.csv'

df = pd.read_csv(filename)

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



# out = adipose_iteration(800, .002)
out = adipose_iteration(800, 1)

# out = adipose_iteration(652, .1)
# out = adipose_iteration(671, .1)

out.to_csv('PHY4003/dose_depth.csv', index=False)


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 13})

fig, ax = plt.subplots()

out['energy_deposited_normalised'] = out['energy_deposited'] / out['energy_deposited'].max()

out.plot(
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