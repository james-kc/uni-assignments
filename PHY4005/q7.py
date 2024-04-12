import numpy as np
from matplotlib import pyplot as plt
from astropy.timeseries import LombScargle


plt.rcParams["font.family"] = "Consolas"
plt.rcParams.update({'font.size': 13})

Observing_times = np.load(f"PHY4005/data/Observing_times.npy")
Boblins_star_BIS = np.load(f"PHY4005/data/Boblins_star_BIS.npy")
Boblins_star_RVs = np.load(f"PHY4005/data/Boblins_star_RVs.npy")
Watsons_star_BIS = np.load(f"PHY4005/data/Watsons_star_BIS.npy")
Watsons_star_RVs = np.load(f"PHY4005/data/Watsons_star_RVs.npy")

Boblins_star_BIS_uncertainties = np.full_like(Boblins_star_BIS, 1)
Boblins_star_RVs_uncertainties = np.full_like(Boblins_star_RVs, 2)
Watsons_star_BIS_uncertainties = np.full_like(Watsons_star_BIS, 0.7)
Watsons_star_RVs_uncertainties = np.full_like(Watsons_star_RVs, 3)


data_files = {
    'Boblins_star_BIS': Boblins_star_BIS,
    'Boblins_star_RVs': Boblins_star_RVs,
    'Watsons_star_BIS': Watsons_star_BIS,
    'Watsons_star_RVs': Watsons_star_RVs
}

data_uncertainties = {
    'Boblins_star_BIS': Boblins_star_BIS_uncertainties,
    'Boblins_star_RVs': Boblins_star_RVs_uncertainties,
    'Watsons_star_BIS': Watsons_star_BIS_uncertainties,
    'Watsons_star_RVs': Watsons_star_RVs_uncertainties
}



for data in data_files:
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(
        Observing_times,
        data_files[data],
        'k'
    )

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
        f'PHY4005/plots/{data}.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()

    # Do not search for periods greater than 50 days since features will be
    # dominated by data gaps and the baseline of the data themselves.
    max_period = 50  # in days
    min_frequency = 1.0 / max_period

    frequency, power = LombScargle(
        Observing_times,
        data_files[data],
        dy=data_uncertainties[data]
    ).autopower(minimum_frequency=min_frequency)

    plt.figure(figsize=(10, 5))
    plt.plot(1 / frequency, power)  # Convert frequency to period
    plt.xlabel('Period')
    plt.ylabel('GLS Power')
    plt.title('Generalized Lomb-Scargle Periodogram')
    plt.grid(True)
    plt.show()