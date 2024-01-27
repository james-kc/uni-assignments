import pandas as pd
from matplotlib import pyplot as plt

filename = 'PHY4003/stopping_power.csv'

df = pd.read_csv(filename)



plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 13})

fig, ax = plt.subplots()
df.plot(
    ax=ax,
    x='kinetic_energy',
    y='stopping_power'
)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(ymin=0)
ax.set_xlim(xmin=0, xmax=30)
ax.get_legend().remove()
ax.set_xlabel(
    'Kinetic Energy (MeV)',
    fontsize=20
)
ax.set_ylabel(
    'Stopping Power (MeV/cm$^2$/g)',
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


plt.show()