import numpy as np
import matplotlib.pyplot as plt

# Constants
epsilon_0 = 8.854187817e-12  # Vacuum permittivity [F/m]
k_B = 1.380649e-23  # Boltzmann constant [J/K]
q_e = 1.602176634e-19  # Elementary charge [C]

# Range of ne and kBT values
ne_values = np.logspace(6, 25, 100)  # Electron number density [m^-3]
kBT_values = np.logspace(-2, 5, 100)  # Temperature [eV]

# Function to calculate Debye length
def debye_length(ne, kBT):
    return np.sqrt((epsilon_0 * kBT) / (ne * q_e**2))

# Calculate Debye length values
lambda_D_values = debye_length(ne_values, kBT_values)

# Plot
plt.figure(figsize=(8, 6))

# Plot lines of constant Debye length
for lambda_D in [1e-3, 1e-2, 1e-1, 1]:
    plt.plot(
        kBT_values,
        (lambda_D/kBT_values)**2,
        linestyle='--',
        color='blue',
        label=f'λ_D={lambda_D}'
    )

# Plot lines of constant density
plt.plot(kBT_values, np.log(kBT_values) - 8, linestyle='--', color='red', label='N_D constant')

# Set log-log scale
plt.xscale('log')
plt.yscale('log')

# Labels and legend
plt.xlabel('kBT (eV)')
plt.ylabel('ne (m^-3)')
plt.title('Log-Log Plot of ne vs kBT with Lines of Constant λ_D and N_D')
plt.legend()
plt.grid(True)

# Show plot
plt.show()
