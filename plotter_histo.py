import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import skew, kurtosis

# --------------------------------------------------
# Load phase space
# --------------------------------------------------

raw = np.loadtxt("Results/ps/phase_space_99.dat")
M = raw.reshape(-1, 6)

x  = M[:,0]   # mm
xp = M[:,1]
y  = M[:,2]   # mm
yp = M[:,3]
t  = M[:,4]   # mm/c
P  = M[:,5]


for name, arr in [('x',x), ('y',y), ('t',t)]:

    print(f"\n{name}")
    print("skewness =", skew(arr))
    print("excess kurtosis =", kurtosis(arr))

# --------------------------------------------------
# Function
# --------------------------------------------------

def gaussian_hist(ax, data, xlabel, bins=80):

    mu = np.mean(data)
    sigma = np.std(data)

    # Histogram normalized to density
    counts, edges, _ = ax.hist(
        data,
        bins=bins,
        density=True,
        color='royalblue',
        alpha=0.7,
        edgecolor='black',
        label='Simulation'
    )

    xfit = np.linspace(edges[0], edges[-1], 1000)

    ax.plot(
        xfit,
        norm.pdf(xfit, mu, sigma),
        'r-',
        lw=3,
        label='Gaussian fit'
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend()

    ax.text(
        0.97,
        0.95,
        rf'$\mu={mu:.3g}$' + '\n' +
        rf'$\sigma={sigma:.3g}$',
        transform=ax.transAxes,
        ha='right',
        va='top',
        bbox=dict(facecolor='white', alpha=0.8)
    )

# --------------------------------------------------
# Plot
# --------------------------------------------------

fig, axs = plt.subplots(1, 3, figsize=(15,5))

gaussian_hist(axs[0], x, 'x [mm]')
gaussian_hist(axs[1], y, 'y [mm]')
gaussian_hist(axs[2], t, 't [mm/c]')

plt.tight_layout()
plt.show()
