import numpy as np
from matplotlib import pyplot as plt
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit



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

data_periods = {
    # 'Boblins_star_BIS': 2.42, # Actual: 2.42 - not correct, coppied from RV
    'Boblins_star_BIS': 27.23, # Actual: 2.42 - not correct, coppied from RV
    'Boblins_star_RVs': 27.23,
    'Watsons_star_BIS': 24.78,
    'Watsons_star_RVs': 24.87
}

def my_uncertainties_plot(
        x,
        y,
        y_uncertainties,
        x_label,
        y_label,
        label,
        filename=None
):
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.errorbar(
        x,
        y,
        yerr=y_uncertainties,
        label=label,
        fmt='o',
        color='blue'
    )

    ax.set_xlabel(
        x_label,
        fontsize=20
    )
    ax.set_ylabel(
        y_label,
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

    ax.legend()

    if filename:
        plt.savefig(
            filename,
            dpi=300,
            bbox_inches='tight'
        )

        plt.show()

    else:
        return ax


def my_plot(x, y, x_label, y_label, filename):

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(x, y, 'k')

    ax.set_xlabel(
        x_label,
        fontsize=20
    )
    ax.set_ylabel(
        y_label,
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
        filename,
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()


# Sine wave
def sinusoidal(t, A, f, phi, c):
    return A * np.sin(2 * np.pi * f * t + phi) + c

for data in data_files:

    ax = my_uncertainties_plot(
        Observing_times,
        data_files[data],
        data_uncertainties[data],
        "Time (days)",
        "Velocity ($ms^{-1}$)",
        'Original Data',
        f'PHY4005/plots/{data}.png'
    )

    # Do not search for periods greater than 50 days since features will be
    # dominated by data gaps and the baseline of the data themselves.
    max_period = 50  # in days
    min_frequency = 1.0 / max_period
    min_period = 1.5  # in days
    max_frequency = 1.0 / min_period

    # Generalised Lomb Scargle periodogram
    frequency, power = LombScargle(
        Observing_times,
        data_files[data],
        dy=data_uncertainties[data]
    ).autopower(
        maximum_frequency=max_frequency,
        minimum_frequency=min_frequency,
        samples_per_peak=1000
    )

    # Find the index of the maximum power
    max_power_index = np.argmax(power)

    # Get the corresponding frequency
    frequency_at_max_power = frequency[max_power_index]

    # Calculate the period from the frequency
    period_at_max_power = 1 / frequency_at_max_power

    print(f"{data}_GLS\t{period_at_max_power} days")

    my_plot(
        1 / frequency,
        power,
        "Period",
        "GLS Power",
        f'PHY4005/plots/{data}_GLS.png'
    )

    # Using periods obtained from periodograms
    period = data_periods[data]

    # Initial guess for the parameters of the sinusoidal function
    A_guess = (max(data_files[data]) - min(data_files[data])) / 2
    phi_guess = 0  # Initial phase guess
    c_guess = np.mean(data_files[data])

    # Use curve_fit to fit the sinusoidal function to the data
    popt, pcov = curve_fit(
        sinusoidal,
        Observing_times,
        data_files[data],
        sigma=data_uncertainties[data],
        p0=[
            A_guess,
            1/period,
            phi_guess,
            c_guess
        ]
    )

    # Plot the original data and the fitted sinusoidal curve
    ax = my_uncertainties_plot(
        Observing_times,
        data_files[data],
        data_uncertainties[data],
        'Time',
        'Radial Velocity',
        'Original Data'
    )
    
    ax.plot(
        Observing_times,
        sinusoidal(Observing_times, *popt),
        label='Fitted Curve',
        color='red'
    )

    ax.legend()

    plt.savefig(
        f'PHY4005/plots/{data}_best_fit.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()

    # Print the parameters of the fitted sinusoidal curve
    print("Amplitude (A):", popt[0])
    print("Frequency (f):", popt[1])
    print("Phase (phi):", popt[2])
    print("Constant (c):", popt[3])
    print()

    # Calculate the phase of each observation
    phase = (Observing_times / period) % 1  # Phase-fold the observations

    # Sort the observations by phase
    sorted_indices = np.argsort(phase)
    phase_sorted = phase[sorted_indices]
    data_sorted = data_files[data][sorted_indices]

    # Plot the phase-folded data points

    ax = my_uncertainties_plot(
        phase_sorted,
        data_sorted,
        data_uncertainties[data][sorted_indices],
        'Phase',
        'Radial Velocity',
        'Phase-folded Data Points'
    )

    # Plot the best-fitting model over the phase-folded data
    ax.plot(phase_sorted, sinusoidal(Observing_times[sorted_indices], *popt), label='Best-fitting Model', color='red')

    ax.legend()

    plt.savefig(
        f'PHY4005/plots/{data}_phase-folded.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()