import numpy as np
from matplotlib import pyplot as plt
import scipy.constants as const


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 13})


def inv_debye_wavelength(kT, debye_length):
    return (const.epsilon_0 * kT) / (debye_length**2 * const.e**2)


def inv_debye_particles(kT, debye_particles):
    return (
        ( (4 * const.pi) / (3 * debye_particles) )**2 *
        ( (const.epsilon_0 * kT) / const.e**2 )**3
    )


steps = 100
kTe_values = np.logspace(-2, 5, 100)  # Temperature [eV]
kT_values = kTe_values * const.e
# n_e = np.logspace(6, 25, steps)  # Electron number density [m^-3]

_, ax = plt.subplots(figsize=(9, 6))

# [1e-3, 1e-2, 1e-1, 1]
for debye_length in [1e-3, 1e-2, 1e-1, 1]:
    ax.plot(
        kTe_values,
        inv_debye_wavelength(kT_values, debye_length),
        linestyle='--',
        # color='blue',
        label=f"$\lambda _D$ = {debye_length}m"
    )

for debye_particles in [1e4, 1e6, 1e8, 1e10]:
    ax.plot(
        kTe_values,
        inv_debye_particles(kT_values, debye_particles),
        linestyle=':',
        # color='blue',
        label=f"$N_D$ = {debye_particles:.1e}"
    )


# Set log-log scale
plt.xscale('log')
plt.yscale('log')

# Labels and legend
plt.title(
    'Log-Log Plot of $n_e$ vs $k_BT$ with Lines of Constant $λ_D$ and $N_D$',
    fontsize=20
)
ax.set_xlabel(
    '$k_BT$ (eV)',
    fontsize=20
)
ax.set_ylabel(
    '$n_e$ ($m^{-3}$)',
    fontsize=20
)
ax.tick_params(
    axis='x',
    labelsize=13,
    which='both',
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

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig(
    f'PHY4008/debye_plasma.png',
    dpi=300,
    bbox_inches='tight'
)

plt.show()