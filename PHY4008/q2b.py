from matplotlib import pyplot as plt
import numpy as np


def dispersion_relation(k, c_s, ion_plas_freq):
    omega_sq = (k**2 * c_s**2) / (1 + ( (k**2 * c_s**2) / (ion_plas_freq**2) ))
    return np.sqrt(omega_sq)


def inverse_exp_decrease(w):
    A = 1
    B = 5
    C = 2.6
    D = 1.2

    return A - C * np.exp(B * (w - D))


k = np.linspace(0, 4, 500)
c_s = 1 # idk, doesn't really matter
ion_plas_freq = 1 # also doesn't matter

_, ax = plt.subplots()

ax.plot(
    k,
    dispersion_relation(k, c_s, ion_plas_freq),
    'k'
)

ax.plot(
    k,
    k,
    'r--'
)

ax.plot(
    k,
    np.full(len(k), ion_plas_freq),
    'r--'
)

ax.set_ylim(0, 1.1)
ax.set_xlim(0, 4)
ax.set_xticks([])
ax.set_yticks([])


ax.set_xlabel(
    'k',
    fontsize=20
)
ax.set_ylabel(
    '$\omega$',
    fontsize=20
)

plt.text(
    0.75, 1.15, '$v_g=c_s$',
    fontsize=20
)

plt.text(
    -.4, .98, '$\Omega_{pi}$',
    fontsize=20
)

plt.savefig(
    f'PHY4008/omega_k_q2b.png',
    dpi=300,
    bbox_inches='tight'
)


w = np.linspace(0, 1.009, 500)

_, ax = plt.subplots()

ax.plot(
    w,
    inverse_exp_decrease(w),
    'k'
)

ax.plot(
    [0, 1.1],
    [0.994, 0.994],
    'r--'
)

ax.plot(
    [1.009, 1.009],
    [0, 1.1],
    'r--'
)

ax.set_ylim(0, 1.1)
ax.set_xlim(0, 1.1)
ax.set_xticks([])
ax.set_yticks([])


ax.set_xlabel(
    '$\omega$',
    fontsize=20
)
ax.set_ylabel(
    '$v_g$',
    fontsize=20
)

plt.text(
    -.1, 0.98, '$c_s$',
    fontsize=20
)
plt.text(
    .97, -0.1, '$\Omega_{pi}$',
    fontsize=20
)

plt.savefig(
    f'PHY4008/vg_omega_q2b.png',
    dpi=300,
    bbox_inches='tight'
)

plt.show()