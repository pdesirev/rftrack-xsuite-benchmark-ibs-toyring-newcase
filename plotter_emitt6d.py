import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import kurtosis
cut = 2.5
def transverse_cut(x,xp,y,yp):
    return (
        (np.abs(x  - np.mean(x))  < cut*np.std(x)) &
        (np.abs(xp - np.mean(xp)) < cut*np.std(xp)) &
        (np.abs(y  - np.mean(y))  < cut*np.std(y)) &
        (np.abs(yp - np.mean(yp)) < cut*np.std(yp))
    )

def longitudinal_cut(t,E):
    return (
        (np.abs(t - np.mean(t)) < cut*np.std(t)) &
        (np.abs(E - np.mean(E)) < cut*np.std(E))
    )

def emit2D(X):
    cov = np.cov(X.T)
    return np.sqrt(np.abs(np.linalg.det(cov)))


def emit4D(X):
    cov = np.cov(X.T)
    return np.abs(np.linalg.det(cov))**0.25


def emit6D(X):
    cov = np.cov(X.T)
    return np.abs(np.linalg.det(cov))**(1/6)

# SETTINGS
# -----------------------------
n_turns = 1000
mass = 0.51099895069  # MeV/c^2 (only if needed for conversions)

emit6d_rf  = np.zeros(n_turns)
emit6d_xs  = np.zeros(n_turns)
emit6d_xs2 = np.zeros(n_turns)
emit4d_rf  = np.zeros(n_turns)
emit4d_xs  = np.zeros(n_turns)
emit4d_xs2 = np.zeros(n_turns)
emitz_rf  = np.zeros(n_turns)
emitz_xs  = np.zeros(n_turns)
emitz_xs2 = np.zeros(n_turns)

# ----------------------------------------------------------
# Shapiro tests
# ----------------------------------------------------------

kurt_x_rf  = np.zeros(n_turns)
kurt_xp_rf = np.zeros(n_turns)
kurt_x_xs  = np.zeros(n_turns)
kurt_xp_xs = np.zeros(n_turns)
kurt_x_xs2  = np.zeros(n_turns)
kurt_xp_xs2 = np.zeros(n_turns)
kurt_t_rf = np.zeros(n_turns)
kurt_E_rf = np.zeros(n_turns)
kurt_t_xs = np.zeros(n_turns)
kurt_E_xs = np.zeros(n_turns)
kurt_t_xs2 = np.zeros(n_turns)
kurt_E_xs2 = np.zeros(n_turns)
def excess_kurtosis(x):
    """
    Excess kurtosis (Fisher definition).

    Returns:
        0  : Gaussian
        >0 : heavier tails than Gaussian
        <0 : lighter tails than Gaussian
    """
    x = np.asarray(x)

    if len(x) < 20:
        return np.nan

    return np.abs(kurtosis(x, fisher=True, bias=False))


for i in range(1,n_turns+1):

    filename_rft = f"Results/ps/phase_space_{i}.dat"
    filename_ps  = f"Results/PS_XSuite_track/ps_{i}.dat"
    filename_xs  = f"Results/PS_XSuite_ana/ps_{i}.dat"

    # ========================================================
    # RF TRACK
    # ========================================================
    M = np.loadtxt(filename_rft).reshape(-1,6)
    x_rft  = M[:,0]
    xp_rft = M[:,1]
    y_rft  = M[:,2]
    yp_rft = M[:,3]
    t_rft = M[:,4] # convert mm/c -> s
    P_rft = M[:,5]
    E_rft = np.sqrt(P_rft**2 + mass**2)
    mask_t_rft = transverse_cut(x_rft,xp_rft,y_rft,yp_rft)
    mask_z_rft = longitudinal_cut(t_rft,E_rft)
    mask_6_rft = mask_t_rft & mask_z_rft

    kurt_x_rf[i-1]  = excess_kurtosis(x_rft[mask_t_rft])
    kurt_xp_rf[i-1] = excess_kurtosis(xp_rft[mask_t_rft])
    kurt_t_rf[i-1]  = excess_kurtosis(t_rft[mask_z_rft])
    kurt_E_rf[i-1]  = excess_kurtosis(E_rft[mask_z_rft])

    X4_rft = np.column_stack([x_rft[mask_t_rft],xp_rft[mask_t_rft],
                          y_rft[mask_t_rft],yp_rft[mask_t_rft]])
    Z2_rft = np.column_stack([t_rft[mask_z_rft],E_rft[mask_z_rft]])
    #X4_rft = np.column_stack([x_rft[mask_6_rft],xp_rft[mask_6_rft],
    #                      y_rft[mask_6_rft],yp_rft[mask_6_rft]])
    #Z2_rft = np.column_stack([t_rft[mask_6_rft],E_rft[mask_6_rft]])
    X6_rft = np.column_stack([x_rft[mask_6_rft],xp_rft[mask_6_rft],
                          y_rft[mask_6_rft],yp_rft[mask_6_rft],
                          t_rft[mask_6_rft],E_rft[mask_6_rft]])

    emit4d_rf[i-1] = emit4D(X4_rft)
    emitz_rf[i-1] = emit2D(Z2_rft)
    emit6d_rf[i-1] = emit6D(X6_rft)
    # ========================================================
    # XSUITE TRACK
    # ========================================================
    M = np.loadtxt(filename_ps)
    x_xt  = M[:,0]
    xp_xt = M[:,1]
    y_xt  = M[:,2]
    yp_xt = M[:,3]
    t_xt  = M[:,4]
    P_xt  = M[:,5]
    E_xt  = np.sqrt(P_xt**2 + mass**2)

    mask_t_xt = transverse_cut(x_xt, xp_xt, y_xt, yp_xt)
    mask_z_xt = longitudinal_cut(t_xt, E_xt)
    mask_6_xt = mask_t_xt & mask_z_xt

    kurt_x_xs[i-1]  = excess_kurtosis(x_xt[mask_t_xt])
    kurt_xp_xs[i-1] = excess_kurtosis(xp_xt[mask_t_xt])
    kurt_t_xs[i-1]  = excess_kurtosis(t_xt[mask_z_xt])
    kurt_E_xs[i-1]  = excess_kurtosis(E_xt[mask_z_xt])

    X4_xt = np.column_stack([
        x_xt[mask_t_xt], xp_xt[mask_t_xt],
        y_xt[mask_t_xt], yp_xt[mask_t_xt]])
    Z2_xt = np.column_stack([
        t_xt[mask_z_xt], E_xt[mask_z_xt]])
    #X4_xt = np.column_stack([
    #    x_xt[mask_6_xt], xp_xt[mask_6_xt],
    #    y_xt[mask_6_xt], yp_xt[mask_6_xt]])
    #Z2_xt = np.column_stack([
    #    t_xt[mask_6_xt], E_xt[mask_6_xt]])
    X6_xt = np.column_stack([
        x_xt[mask_6_xt], xp_xt[mask_6_xt],
        y_xt[mask_6_xt], yp_xt[mask_6_xt],
        t_xt[mask_6_xt], E_xt[mask_6_xt]])

    emit4d_xs[i-1] = emit4D(X4_xt)
    emitz_xs[i-1] = emit2D(Z2_xt)
    emit6d_xs[i-1] = emit6D(X6_xt)

    # ========================================================
    # XSUITE ANALYTICAL
    # ========================================================
    M = np.loadtxt(filename_xs)
    x_xa  = M[:,0]
    xp_xa = M[:,1]
    y_xa  = M[:,2]
    yp_xa = M[:,3]
    t_xa  = M[:,4]
    P_xa  = M[:,5]
    E_xa  = np.sqrt(P_xa**2 + mass**2)

    mask_t_xa = transverse_cut(x_xa, xp_xa, y_xa, yp_xa)
    mask_z_xa = longitudinal_cut(t_xa, E_xa)
    mask_6_xa = mask_t_xa & mask_z_xa

    kurt_x_xs2[i-1]  = excess_kurtosis(x_xa[mask_t_xa])
    kurt_xp_xs2[i-1] = excess_kurtosis(xp_xa[mask_t_xa])
    kurt_t_xs2[i-1]  = excess_kurtosis(t_xa[mask_z_xa])
    kurt_E_xs2[i-1]  = excess_kurtosis(E_xa[mask_z_xa])

    X4_xa = np.column_stack([
        x_xa[mask_t_xa], xp_xa[mask_t_xa],
        y_xa[mask_t_xa], yp_xa[mask_t_xa]])
    Z2_xa = np.column_stack([
        t_xa[mask_z_xa], E_xa[mask_z_xa]])
    #X4_xa = np.column_stack([
    #    x_xa[mask_6_xa], xp_xa[mask_6_xa],
    #    y_xa[mask_6_xa], yp_xa[mask_6_xa]])
    #Z2_xa = np.column_stack([
    #    t_xa[mask_6_xa], E_xa[mask_6_xa]])
    X6_xa = np.column_stack([
        x_xa[mask_6_xa], xp_xa[mask_6_xa],
        y_xa[mask_6_xa], yp_xa[mask_6_xa],
        t_xa[mask_6_xa], E_xa[mask_6_xa]])

    emit4d_xs2[i-1] = emit4D(X4_xa)
    emitz_xs2[i-1] = emit2D(Z2_xa)
    emit6d_xs2[i-1] = emit6D(X6_xa)


######################################
## WINDOWS
######################################
window_trans = 41
window_long = 11

def moving_average(x, window=20):
    return np.convolve(x, np.ones(window)/window, mode='valid')

kurt_x_xs_ma   = moving_average(kurt_x_xs, window_trans)
kurt_xp_xs_ma  = moving_average(kurt_xp_xs, window_trans)
kurt_t_xs_ma = moving_average(kurt_t_xs, window_long)
kurt_E_xs_ma = moving_average(kurt_E_xs, window_long)

kurt_x_xs2_ma  = moving_average(kurt_x_xs2, window_trans)
kurt_xp_xs2_ma = moving_average(kurt_xp_xs2, window_trans)
kurt_t_xs2_ma = moving_average(kurt_t_xs2, window_long)
kurt_E_xs2_ma = moving_average(kurt_E_xs2, window_long)

kurt_x_rf_ma   = moving_average(kurt_x_rf, window_trans)
kurt_xp_rf_ma  = moving_average(kurt_xp_rf, window_trans)
kurt_t_rf_ma = moving_average(kurt_t_rf, window_long)
kurt_E_rf_ma = moving_average(kurt_E_rf, window_long)

emit4d_rf_ma  = moving_average(emit4d_rf, window_trans)
emit4d_xs_ma  = moving_average(emit4d_xs, window_trans)
emit4d_xs2_ma = moving_average(emit4d_xs2, window_trans)

emitz_rf_ma  = moving_average(emitz_rf, window_long)
emitz_xs_ma  = moving_average(emitz_xs, window_long)
emitz_xs2_ma = moving_average(emitz_xs2, window_long)

# ============================================================
# PLOT COMPARISON
# ============================================================

turns = np.arange(1, n_turns + 1)

# -----------------------------
# 6D emittance
# -----------------------------
plt.figure(figsize=(8,5))

plt.plot(turns, emit6d_rf/emit6d_rf[0], 'r', lw=2, label="RF-Track")
plt.plot(turns, emit6d_xs/emit6d_xs[0], 'b', lw=2, label="Zampetakis et al.")
plt.plot(turns, emit6d_xs2/emit6d_xs2[0], 'k', lw=2, label="Nagaitsev")
plt.legend(fontsize=14, loc='upper left')
plt.grid(True)
plt.xlabel("Number of turns", fontsize=16)
plt.ylabel(r"$\epsilon_{6D}/\epsilon_{6D,0}$", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("6d_emittance.png", dpi=300, bbox_inches="tight")
plt.savefig("6d_emittance.pdf", bbox_inches="tight")
plt.close()


turns_ma = turns[:len(emit4d_rf_ma)] + window_trans//2
turns_maz = turns[:len(emitz_rf_ma)] + window_long//2

plt.figure(figsize=(8,5))
plt.plot(turns_ma, emit4d_rf_ma, 'r', lw=2, label="RF-Track")
plt.plot(turns_ma, emit4d_xs_ma, 'b', lw=2, label="Zampetakis et al.")
plt.plot(turns_ma, emit4d_xs2_ma, 'k', lw=2, label="Nagaitsev")
plt.legend(fontsize=14, loc='upper left')
plt.grid(True)
plt.xlabel("Number of turns", fontsize=16)
plt.ylabel(r"$\epsilon_{4D}/\epsilon_{4D,0}$", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("4d_transverse_emittance_ma20.png", dpi=300, bbox_inches="tight")
plt.savefig("4d_transverse_emittance_ma20.pdf", bbox_inches="tight")
plt.close()


# -----------------------------
# Longitudinal emittance
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(turns_maz, emitz_rf_ma, 'r', lw=2, label="RF-Track")
plt.plot(turns_maz, emitz_xs_ma, 'b', lw=2, label="Zampetakis et al.")
plt.plot(turns_maz, emitz_xs2_ma, 'k', lw=2, label="Nagaitsev")
plt.legend(fontsize=14, loc='upper left')
plt.grid(True)
plt.xlabel("Number of turns", fontsize=16)
plt.ylabel(r"$\epsilon_{z}/\epsilon_{z,0}$", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("longitudinal_emittance.png", dpi=300, bbox_inches="tight")
plt.savefig("longitudinal_emittance.pdf", bbox_inches="tight")
plt.close()

########## TEST ##########
plt.figure(figsize=(8,5))
plt.plot(turns_ma, kurt_x_rf_ma,   'r-',  lw=2, label='RF-Track x')
plt.plot(turns_ma, kurt_xp_rf_ma,  'r--', lw=2, label="RF-Track x'")
plt.plot(turns_ma, kurt_x_xs_ma,   'b-',  lw=2, label='Xsuite x')
plt.plot(turns_ma, kurt_xp_xs_ma,  'b--', lw=2, label="Xsuite x'")
plt.plot(turns_ma, kurt_x_xs2_ma,  'k-',  lw=2, label='Nagaitsev x')
plt.plot(turns_ma, kurt_xp_xs2_ma, 'k--', lw=2, label="Nagaitsev x'")
plt.grid(True)
plt.xlabel("Turn", fontsize=16)
plt.ylabel("Excess kurtosis", fontsize=16)
plt.legend(fontsize=11, ncol=2)
plt.tight_layout()
plt.savefig("kurtosis_transverse.png", dpi=300)
plt.savefig("kurtosis_transverse.pdf")
plt.close()

plt.figure(figsize=(8,5))
plt.plot(turns_maz, kurt_t_rf_ma, 'r-', lw=2, label='RF-Track t')
plt.plot(turns_maz, kurt_E_rf_ma, 'r--', lw=2, label='RF-Track E')
plt.plot(turns_maz, kurt_t_xs_ma, 'b-', lw=2, label='Xsuite t')
plt.plot(turns_maz, kurt_E_xs_ma, 'b--', lw=2, label='Xsuite E')
plt.plot(turns_maz, kurt_t_xs2_ma, 'k-', lw=2, label='Nagaitsev t')
plt.plot(turns_maz, kurt_E_xs2_ma, 'k--', lw=2, label='Nagaitsev E')
plt.grid(True)
plt.xlabel("Turn", fontsize=16)
plt.ylabel("Excess kurtosis", fontsize=16)
plt.legend(fontsize=11, ncol=2)
plt.tight_layout()
plt.savefig("kurtosis_longitudinal.png", dpi=300)
plt.savefig("kurtosis_longitudinal.pdf")
plt.close()