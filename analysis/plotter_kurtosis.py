import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


# ============================================================
# SETTINGS
# ============================================================

n_turns_rf = 4000
n_files_xs = 160

# Xsuite phase-space files are saved every 25 turns
xs_turn_spacing = 25

# Moving-average windows
window_trans = 41   # RF-Track transverse
window_long  = 11   # RF-Track longitudinal
window_xs    = 3    # Xsuite: very light smoothing

# Expected two-sided Gaussian population outside +/- 3 sigma
gaussian_tail = 0.0026998


# ============================================================
# FUNCTIONS
# ============================================================

def robust_tail_fraction(x, threshold=3):
    """
    Fraction of particles farther than 'threshold' robust
    standard deviations from the median.

    sigma_rob = 1.4826 * MAD

    where

        MAD = median(|x - median(x)|)

    For a Gaussian distribution and threshold = 3,
    the expected two-sided tail population is ~0.27%.
    """

    x = np.asarray(x)
    x = x[np.isfinite(x)]

    med = np.median(x)
    mad = np.median(np.abs(x - med))

    robust_sigma = 1.4826 * mad

    if robust_sigma == 0:
        return np.nan

    return np.mean(
        np.abs(x - med) > threshold * robust_sigma
    )


def moving_average(x, window):
    """
    Moving average using a fixed window.
    """

    x = np.asarray(x)

    if len(x) < window:
        return x

    return np.convolve(
        x,
        np.ones(window) / window,
        mode='valid'
    )


# ============================================================
# ARRAYS
# ============================================================

# RF-Track
tail_x_rf  = []
tail_xp_rf = []
tail_t_rf  = []
tail_p_rf  = []

# Xsuite - Zampetakis et al.
tail_x_xs  = []
tail_xp_xs = []
tail_t_xs  = []
tail_p_xs  = []

# Xsuite - Nagaitsev
tail_x_nag  = []
tail_xp_nag = []
tail_t_nag  = []
tail_p_nag  = []


# ============================================================
# COMMON INITIAL DISTRIBUTION -- TURN 0
# ============================================================
#
# initial_bunch_25pC.txt follows the RF-Track input convention:
#
# column 0 : x
# column 1 : x'
# column 2 : y
# column 3 : y'
# column 4 : t
# column 5 : p [MeV/c]
# column 6 : mass
# column 7 : charge
# column 8 : weight
#
# The same physical particle distribution is used as the
# initial condition for RF-Track and both Xsuite calculations.
# ============================================================

M0 = np.loadtxt("initial_bunch_25pC.txt")

x0  = M0[:, 0]
xp0 = M0[:, 1]
t0  = M0[:, 4]
p0  = M0[:, 5]


# Initial tail populations
tail_x0  = robust_tail_fraction(x0)
tail_xp0 = robust_tail_fraction(xp0)
tail_t0  = robust_tail_fraction(t0)
tail_p0  = robust_tail_fraction(p0)


print("\n============================================")
print("COMMON INITIAL DISTRIBUTION -- TURN 0")
print("============================================")
print(f"x  : {100 * tail_x0:.3f} %")
print(f"x' : {100 * tail_xp0:.3f} %")
print(f"t  : {100 * tail_t0:.3f} %")
print(f"p  : {100 * tail_p0:.3f} %")
print("============================================\n")


# ============================================================
# ADD TURN 0 TO ALL METHODS
# ============================================================

# RF-Track
tail_x_rf.append(tail_x0)
tail_xp_rf.append(tail_xp0)
tail_t_rf.append(tail_t0)
tail_p_rf.append(tail_p0)

# Xsuite - Zampetakis
tail_x_xs.append(tail_x0)
tail_xp_xs.append(tail_xp0)
tail_t_xs.append(tail_t0)
tail_p_xs.append(tail_p0)

# Xsuite - Nagaitsev
tail_x_nag.append(tail_x0)
tail_xp_nag.append(tail_xp0)
tail_t_nag.append(tail_t0)
tail_p_nag.append(tail_p0)


# ============================================================
# RF-TRACK
# ============================================================
#
# RF-Track convention:
#
# turn 0 -> initial_bunch_25pC.txt
# turn 1 -> phase_space_1.dat
# turn 2 -> phase_space_2.dat
# ...
# ============================================================

for i in range(1, n_turns_rf + 1):

    filename_rft = (
        f"../Results/ps/"
        f"phase_space_{i}.dat"
    )

    M = np.loadtxt(filename_rft).reshape(-1, 6)

    x_rft  = M[:, 0]
    xp_rft = M[:, 1]
    t_rft  = M[:, 4]
    p_rft  = M[:, 5]      # MeV/c

    tail_x_rf.append(
        robust_tail_fraction(x_rft)
    )

    tail_xp_rf.append(
        robust_tail_fraction(xp_rft)
    )

    tail_t_rf.append(
        robust_tail_fraction(t_rft)
    )

    tail_p_rf.append(
        robust_tail_fraction(p_rft)
    )


# ============================================================
# XSUITE
# ============================================================
#
# Xsuite convention used here:
#
# turn 0  -> initial_bunch_25pC.txt
# ps_1    -> turn 25
# ps_2    -> turn 50
# ps_3    -> turn 75
# ...
# ============================================================

for i in range(1, n_files_xs + 1):

    # --------------------------------------------------------
    # ZAMPETAKIS ET AL.
    # --------------------------------------------------------

    filename_ps = (
        f"../Results/PS_XSuite_track/"
        f"ps_{i}.dat"
    )

    M = np.loadtxt(filename_ps)

    x_xt  = M[:, 0]
    xp_xt = M[:, 1]
    t_xt  = M[:, 4]
    p_xt  = M[:, 5]      # MeV/c

    tail_x_xs.append(
        robust_tail_fraction(x_xt)
    )

    tail_xp_xs.append(
        robust_tail_fraction(xp_xt)
    )

    tail_t_xs.append(
        robust_tail_fraction(t_xt)
    )

    tail_p_xs.append(
        robust_tail_fraction(p_xt)
    )


    # --------------------------------------------------------
    # NAGAITSEV
    # --------------------------------------------------------

    filename_nag = (
        f"../Results/PS_XSuite_ana/"
        f"ps_{i}.dat"
    )

    M = np.loadtxt(filename_nag)

    x_nag  = M[:, 0]
    xp_nag = M[:, 1]
    t_nag  = M[:, 4]
    p_nag  = M[:, 5]      # MeV/c

    tail_x_nag.append(
        robust_tail_fraction(x_nag)
    )

    tail_xp_nag.append(
        robust_tail_fraction(xp_nag)
    )

    tail_t_nag.append(
        robust_tail_fraction(t_nag)
    )

    tail_p_nag.append(
        robust_tail_fraction(p_nag)
    )


# ============================================================
# CONVERT TO NUMPY ARRAYS
# ============================================================

tail_x_rf  = np.asarray(tail_x_rf)
tail_xp_rf = np.asarray(tail_xp_rf)
tail_t_rf  = np.asarray(tail_t_rf)
tail_p_rf  = np.asarray(tail_p_rf)

tail_x_xs  = np.asarray(tail_x_xs)
tail_xp_xs = np.asarray(tail_xp_xs)
tail_t_xs  = np.asarray(tail_t_xs)
tail_p_xs  = np.asarray(tail_p_xs)

tail_x_nag  = np.asarray(tail_x_nag)
tail_xp_nag = np.asarray(tail_xp_nag)
tail_t_nag  = np.asarray(tail_t_nag)
tail_p_nag  = np.asarray(tail_p_nag)


# ============================================================
# TURN ARRAYS
# ============================================================

# RF-Track:
# 0, 1, 2, ..., 4000
turns_rf = np.arange(
    0,
    n_turns_rf + 1
)


# Xsuite:
# 0, 25, 50, ..., 4000
turns_xs = np.arange(
    0,
    (n_files_xs + 1) * xs_turn_spacing,
    xs_turn_spacing
)


# ============================================================
# CHECK ARRAY LENGTHS
# ============================================================

print("Array lengths:")
print(
    "RF-Track:",
    len(turns_rf),
    len(tail_x_rf)
)

print(
    "Xsuite:",
    len(turns_xs),
    len(tail_x_xs)
)

print(
    "Nagaitsev:",
    len(turns_xs),
    len(tail_x_nag)
)

print()


# ============================================================
# PRINT EARLY RF-TRACK EVOLUTION
# ============================================================

print("RF-TRACK TAIL POPULATION")
print("================================================")
print("Turn       x       x'       t       p")
print("------------------------------------------------")

for turn in [
    0,
    1,
    5,
    10,
    20,
    25,
    50,
    100,
    250,
    500
]:

    print(
        f"{turn:4d}   "
        f"{100 * tail_x_rf[turn]:6.3f}%  "
        f"{100 * tail_xp_rf[turn]:6.3f}%  "
        f"{100 * tail_t_rf[turn]:6.3f}%  "
        f"{100 * tail_p_rf[turn]:6.3f}%"
    )

print("================================================\n")


# ============================================================
# RF-TRACK MOVING AVERAGES
# ============================================================

# ------------------------------------------------------------
# Transverse -- 41-turn window
# ------------------------------------------------------------

tail_x_rf_ma = moving_average(
    tail_x_rf,
    window_trans
)

tail_xp_rf_ma = moving_average(
    tail_xp_rf,
    window_trans
)

# Applying the same moving average to the turn coordinate
# places each value at the centre of its averaging interval.
turns_rf_trans_ma = moving_average(
    turns_rf,
    window_trans
)


# ------------------------------------------------------------
# Longitudinal -- 11-turn window
# ------------------------------------------------------------

tail_t_rf_ma = moving_average(
    tail_t_rf,
    window_long
)

tail_p_rf_ma = moving_average(
    tail_p_rf,
    window_long
)

turns_rf_long_ma = moving_average(
    turns_rf,
    window_long
)


# ============================================================
# XSUITE MOVING AVERAGES
# ============================================================
#
# Very light smoothing:
#
# window_xs = 2
#
# Since Xsuite is sampled every 25 turns, each plotted value
# averages two consecutive stored distributions.
# ============================================================

# ------------------------------------------------------------
# Zampetakis et al.
# ------------------------------------------------------------

tail_x_xs_ma = moving_average(
    tail_x_xs,
    window_xs
)

tail_xp_xs_ma = moving_average(
    tail_xp_xs,
    window_xs
)

tail_t_xs_ma = moving_average(
    tail_t_xs,
    window_xs
)

tail_p_xs_ma = moving_average(
    tail_p_xs,
    window_xs
)


# ------------------------------------------------------------
# Nagaitsev
# ------------------------------------------------------------

tail_x_nag_ma = moving_average(
    tail_x_nag,
    window_xs
)

tail_xp_nag_ma = moving_average(
    tail_xp_nag,
    window_xs
)

tail_t_nag_ma = moving_average(
    tail_t_nag,
    window_xs
)

tail_p_nag_ma = moving_average(
    tail_p_nag,
    window_xs
)


# ------------------------------------------------------------
# Xsuite turn coordinates
# ------------------------------------------------------------

turns_xs_ma = moving_average(
    turns_xs,
    window_xs
)


# ============================================================
# TRANSVERSE TAIL POPULATION
# ============================================================

plt.figure(
    figsize=(8, 5)
)


# RF-Track
plt.plot(
    turns_rf_trans_ma,
    tail_x_rf_ma,
    'r-',
    linewidth=2,
    label='RF-Track - x'
)

plt.plot(
    turns_rf_trans_ma,
    tail_xp_rf_ma,
    'r--',
    linewidth=2,
    label="RF-Track - x'"
)


# Xsuite - Zampetakis
plt.plot(
    turns_xs_ma,
    tail_x_xs_ma,
    'b-',
    linewidth=2,
    label='Zampetakis et al. (Xsuite) - x'
)

plt.plot(
    turns_xs_ma,
    tail_xp_xs_ma,
    'b--',
    linewidth=2,
    label="Zampetakis et al. (Xsuite) - x'"
)


# Xsuite - Nagaitsev
plt.plot(
    turns_xs_ma,
    tail_x_nag_ma,
    'k-',
    linewidth=2,
    label='Nagaitsev (Xsuite) - x'
)

plt.plot(
    turns_xs_ma,
    tail_xp_nag_ma,
    'k--',
    linewidth=2,
    label="Nagaitsev (Xsuite) - x'"
)


# Gaussian reference
plt.axhline(
    gaussian_tail,
    color='0.45',
    linestyle=':',
    linewidth=1.5,
    label='Gaussian reference'
)


# Axis formatting
plt.gca().yaxis.set_major_formatter(
    PercentFormatter(
        xmax=1.0,
        decimals=2
    )
)

# Natural zero baseline; upper limit = 0.50 %
plt.ylim(
    0.0,
    0.01
)

plt.xlim(
    0,
    n_turns_rf
)

plt.grid(True)

plt.xlabel(
    "Turn",
    fontsize=16
)

plt.ylabel(
    "Tail population",
    fontsize=16
)

plt.legend(
    fontsize=10,
    ncol=2
)

plt.tight_layout()


plt.savefig(
    "tail_population_transverse.png",
    dpi=300
)

plt.savefig(
    "tail_population_transverse.pdf"
)

plt.close()


# ============================================================
# LONGITUDINAL TAIL POPULATION
# ============================================================

plt.figure(
    figsize=(8, 5)
)


# RF-Track
plt.plot(
    turns_rf_long_ma,
    tail_t_rf_ma,
    'r-',
    linewidth=2,
    label='RF-Track - t'
)

plt.plot(
    turns_rf_long_ma,
    tail_p_rf_ma,
    'r--',
    linewidth=2,
    label='RF-Track - p'
)


# Xsuite - Zampetakis
plt.plot(
    turns_xs_ma,
    tail_t_xs_ma,
    'b-',
    linewidth=2,
    label='Zampetakis et al. (Xsuite) - t'
)

plt.plot(
    turns_xs_ma,
    tail_p_xs_ma,
    'b--',
    linewidth=2,
    label='Zampetakis et al. (Xsuite) - p'
)


# Xsuite - Nagaitsev
plt.plot(
    turns_xs_ma,
    tail_t_nag_ma,
    'k-',
    linewidth=2,
    label='Nagaitsev (Xsuite) - t'
)

plt.plot(
    turns_xs_ma,
    tail_p_nag_ma,
    'k--',
    linewidth=2,
    label='Nagaitsev (Xsuite) - p'
)


# Gaussian reference
plt.axhline(
    gaussian_tail,
    color='0.45',
    linestyle=':',
    linewidth=1.5,
    label='Gaussian reference'
)


# Axis formatting
plt.gca().yaxis.set_major_formatter(
    PercentFormatter(
        xmax=1.0,
        decimals=2
    )
)

# Natural zero baseline; upper limit = 1.80 %
plt.ylim(
    0.0,
    0.018
)

plt.xlim(
    0,
    n_turns_rf
)

plt.grid(True)

plt.xlabel(
    "Turn",
    fontsize=16
)

plt.ylabel(
    "Tail population",
    fontsize=16
)

plt.legend(
    fontsize=10,
    ncol=2
)

plt.tight_layout()


plt.savefig(
    "tail_population_longitudinal.png",
    dpi=300
)

plt.savefig(
    "tail_population_longitudinal.pdf"
)

plt.close()