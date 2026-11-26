import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import re
from scipy.stats import kurtosis
from sklearn.covariance import MinCovDet
cut_trans = 3
cut_long = 3

def mahalanobis_cut(X, cut):
    mean = np.mean(X, axis=0)
    cov = np.cov(X.T)
    cov += np.eye(cov.shape[0])*1e-12
    cov_inv = np.linalg.inv(cov)
    d = X - mean
    R2 = np.einsum(
        'ij,jk,ik->i',
        d,
        cov_inv,
        d
    )
    return R2 < cut**2

def transverse_cut_mahalanobis(x,xp,y,yp,cut):

    X = np.column_stack([
        x,
        xp,
        y,
        yp
    ])
    return mahalanobis_cut(X,cut)

def plane_cut_mahalanobis(a, ap, cut):

    X = np.column_stack([a, ap])

    return mahalanobis_cut(X, cut)

def longitudinal_cut(t,E, cut):
    t0 = np.mean(t)
    E0 = np.mean(E)

    st = 1.4826 * np.median(np.abs(t - t0))
    sE = 1.4826 * np.median(np.abs(E - E0))

    R2 = ((t - t0)/st)**2 + ((E - E0)/sE)**2

    return R2 < cut**2

def emit2D(X):
    cov = np.cov(X.T)
    return np.sqrt(np.abs(np.linalg.det(cov)))

def emit4D(X):
    cov = np.cov(X.T)
    return np.abs(np.linalg.det(cov))**0.25


def emit6D(X):
    cov = np.cov(X.T)
    return np.abs(np.linalg.det(cov))**(1/6)

mass = 0.51099895069  # MeV/c^2 (only if needed for conversions)

emit6d_rf  = []
emit6d_xs  = []
emit6d_xs2 = []
emit4d_rf  = []
emit4d_xs  = []
emit4d_xs2 = []
emitz_rf  = []
emitz_xs  = []
emitz_xs2 = []

emitx_rf  = []
emity_rf  = []
emitx_xs  = []
emity_xs  = []
emitx_xs2 = []
emity_xs2 = []

# ----------------------------------------------------------
# Shapiro tests
# ----------------------------------------------------------

kurt_x_rf  = []
kurt_xp_rf = []
kurt_x_xs  = []
kurt_xp_xs = []
kurt_x_xs2  = []
kurt_xp_xs2 = []
kurt_t_rf = []
kurt_E_rf = []
kurt_t_xs = []
kurt_E_xs = []
kurt_t_xs2 = []
kurt_E_xs2 = []
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


# Discover RF-Track files dynamically; retain their actual turn indices.
rf_dir = Path("../Results/ps")
max_turns_rf = 4000

rf_files = []
for path in rf_dir.glob("phase_space_*.dat"):
    match = re.fullmatch(r"phase_space_(\d+)\.dat", path.name)
    if match:
        turn = int(match.group(1))
        if turn <= max_turns_rf:
            rf_files.append((turn, path))

rf_files.sort(key=lambda item: item[0])

if not rf_files:
    raise FileNotFoundError(
        f"No phase_space_*.dat files found in {rf_dir} up to turn {max_turns_rf}"
    )

# Skip zero-byte files left by incomplete or older runs.
nonempty_rf_files = [
    (turn, path) for turn, path in rf_files
    if path.stat().st_size > 0
]

print(
    f"RF-Track: found {len(rf_files)} files up to turn {max_turns_rf}; "
    f"{len(rf_files) - len(nonempty_rf_files)} are empty"
)

rf_files = nonempty_rf_files

if not rf_files:
    raise ValueError("No non-empty RF-Track phase-space files found")

rf_turns_processed = []

n_turns_xs = 160       # Xsuite tracking available turns (unchanged)

for i, filename_rft in rf_files:

    # ========================================================
    # RF TRACK
    # ========================================================
    try:
        M = np.loadtxt(filename_rft)
        if M.size == 0 or M.size % 6 != 0:
            print(f"Skipping empty or malformed file: {filename_rft}")
            continue
        M = M.reshape(-1, 6)
        if len(M) < 10 or not np.all(np.isfinite(M)):
            print(f"Skipping insufficient/non-finite data: {filename_rft}")
            continue
    except (OSError, ValueError) as exc:
        print(f"Skipping unreadable file {filename_rft}: {exc}")
        continue
    rf_turns_processed.append(i)
    x_rft  = M[:,0]
    xp_rft = M[:,1]
    y_rft  = M[:,2]
    yp_rft = M[:,3]
    t_rft = M[:,4] * 3.335641e-09 # convert mm/c -> s
    P_rft = M[:,5]  # MeV/c
    E_rft = np.sqrt(P_rft**2 + mass**2) * 1e6 

    mask_t_rft = transverse_cut_mahalanobis(x_rft, xp_rft, y_rft, yp_rft, cut_trans)
    mask_z_rft = longitudinal_cut(t_rft,E_rft, cut_long)
    mask_6_rft = mask_t_rft & mask_z_rft

    #mask_x_rft = plane_cut_mahalanobis(x_rft,xp_rft,cut_trans)
    #mask_y_rft = plane_cut_mahalanobis(y_rft,yp_rft,cut_trans)

    #kurt_x_rf.append(excess_kurtosis(x_rft[mask_t_rft]))
    #kurt_xp_rf.append(excess_kurtosis(xp_rft[mask_t_rft]))
    #kurt_t_rf.append(excess_kurtosis(t_rft[mask_z_rft]))
    #kurt_E_rf.append(excess_kurtosis(E_rft[mask_z_rft]))

    kurt_x_rf.append(excess_kurtosis(x_rft))
    kurt_xp_rf.append(excess_kurtosis(xp_rft))
    kurt_t_rf.append(excess_kurtosis(t_rft))
    kurt_E_rf.append(excess_kurtosis(E_rft))

    X4_rft = np.column_stack([x_rft[mask_t_rft],xp_rft[mask_t_rft],
                          y_rft[mask_t_rft],yp_rft[mask_t_rft]])
    Z2_rft = np.column_stack([t_rft[mask_z_rft],E_rft[mask_z_rft]])
    #X4_rft = np.column_stack([x_rft[mask_6_rft],xp_rft[mask_6_rft],
    #                      y_rft[mask_6_rft],yp_rft[mask_6_rft]])
    #Z2_rft = np.column_stack([t_rft[mask_6_rft],E_rft[mask_6_rft]])
    X6_rft = np.column_stack([x_rft[mask_6_rft],xp_rft[mask_6_rft],
                          y_rft[mask_6_rft],yp_rft[mask_6_rft],
                          t_rft[mask_6_rft],E_rft[mask_6_rft]])
    X2x_rft = np.column_stack([
    x_rft[mask_t_rft],
    xp_rft[mask_t_rft]])

    X2y_rft = np.column_stack([
        y_rft[mask_t_rft],
        yp_rft[mask_t_rft]])


    emit4d_rf.append(emit4D(X4_rft))
    emitz_rf.append(emit2D(Z2_rft))
    emit6d_rf.append(emit6D(X6_rft))
    emitx_rf.append(emit2D(X2x_rft))
    emity_rf.append(emit2D(X2y_rft))

rf_turns = np.asarray(rf_turns_processed, dtype=int)
if len(rf_turns) == 0:
    raise ValueError("No usable RF-Track phase-space files found")
print(f"RF-Track: processed {len(rf_turns)} usable files, turns {rf_turns[0]} to {rf_turns[-1]}")

for i in range(1,n_turns_xs+1):
    filename_ps  = f"../Results/PS_XSuite_track/ps_{i}.dat"
    filename_xs  = f"../Results/PS_XSuite_ana/ps_{i}.dat"
    # ========================================================
    # XSUITE TRACK
    # ========================================================
    M = np.loadtxt(filename_ps)
    x_xt  = M[:,0]
    xp_xt = M[:,1]
    y_xt  = M[:,2]
    yp_xt = M[:,3]
    t_xt  = M[:,4] * 3.335641e-09
    P_xt  = M[:,5]
    E_xt  = np.sqrt(P_xt**2 + mass**2) * 1e6

    mask_t_xt = transverse_cut_mahalanobis(x_xt,xp_xt,y_xt,yp_xt,cut_trans)
    mask_z_xt = longitudinal_cut(t_xt,E_xt,cut_long)
    #mask_x_xt = plane_cut_mahalanobis(x_xt, xp_xt, cut_trans)
    #mask_y_xt = plane_cut_mahalanobis(y_xt,yp_xt, cut_trans)
    mask_6_xt = mask_t_xt & mask_z_xt

    #kurt_x_xs.append(excess_kurtosis(x_xt[mask_t_xt]))
    #kurt_xp_xs.append(excess_kurtosis(xp_xt[mask_t_xt]))
    #kurt_t_xs.append(excess_kurtosis(t_xt[mask_z_xt]))
    #kurt_E_xs.append(excess_kurtosis(E_xt[mask_z_xt]))
    kurt_x_xs.append(excess_kurtosis(x_xt))
    kurt_xp_xs.append(excess_kurtosis(xp_xt))
    kurt_t_xs.append(excess_kurtosis(t_xt))
    kurt_E_xs.append(excess_kurtosis(E_xt))

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
    
    X2x_xt = np.column_stack([
        x_xt[mask_t_xt],
        xp_xt[mask_t_xt]])

    X2y_xt = np.column_stack([
        y_xt[mask_t_xt],
        yp_xt[mask_t_xt]])

    emit4d_xs.append(emit4D(X4_xt))
    emitz_xs.append(emit2D(Z2_xt))
    emit6d_xs.append(emit6D(X6_xt))
    emitx_xs.append(emit2D(X2x_xt))
    emity_xs.append(emit2D(X2y_xt))
    # ========================================================
    # XSUITE ANALYTICAL
    # ========================================================
    M = np.loadtxt(filename_xs)
    x_xa  = M[:,0]
    xp_xa = M[:,1]
    y_xa  = M[:,2]
    yp_xa = M[:,3]
    t_xa  = M[:,4] * 3.335641e-09
    P_xa  = M[:,5]
    E_xa  = np.sqrt(P_xa**2 + mass**2) * 1e6

    mask_t_xa = transverse_cut_mahalanobis(
        x_xa,
        xp_xa,
        y_xa,
        yp_xa,
        cut_trans
    )

    mask_z_xa = longitudinal_cut(
        t_xa,
        E_xa,
        cut_long
    )

    mask_6_xa = mask_t_xa & mask_z_xa
    #kurt_x_xs2.append(excess_kurtosis(x_xa[mask_t_xa]))
    #kurt_xp_xs2.append(excess_kurtosis(xp_xa[mask_t_xa]))
    #kurt_t_xs2.append(excess_kurtosis(t_xa[mask_z_xa]))
    #kurt_E_xs2.append(excess_kurtosis(E_xa[mask_z_xa]))
    kurt_x_xs2.append(excess_kurtosis(x_xa))
    kurt_xp_xs2.append(excess_kurtosis(xp_xa))
    kurt_t_xs2.append(excess_kurtosis(t_xa))
    kurt_E_xs2.append(excess_kurtosis(E_xa))


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

    X2x_xa = np.column_stack([
        x_xa[mask_t_xa], xp_xa[mask_t_xa]])

    X2y_xa = np.column_stack([y_xa[mask_t_xa],yp_xa[mask_t_xa]])

    emit4d_xs2.append(emit4D(X4_xa))
    emitz_xs2.append(emit2D(Z2_xa))
    emit6d_xs2.append(emit6D(X6_xa))
    emitx_xs2.append(emit2D(X2x_xa))
    emity_xs2.append(emit2D(X2y_xa))


######################################
## WINDOWS
######################################

emit6d_rf  = np.array(emit6d_rf)
emit6d_xs  = np.array(emit6d_xs)
emit6d_xs2 = np.array(emit6d_xs2)

emit4d_rf  = np.array(emit4d_rf)
emit4d_xs  = np.array(emit4d_xs)
emit4d_xs2 = np.array(emit4d_xs2)

emitz_rf  = np.array(emitz_rf)
emitz_xs  = np.array(emitz_xs)
emitz_xs2 = np.array(emitz_xs2)

emitx_rf  = np.array(emitx_rf)
emity_rf  = np.array(emity_rf)
emitx_xs  = np.array(emitx_xs)
emity_xs  = np.array(emity_xs)
emitx_xs2 = np.array(emitx_xs2)
emity_xs2 = np.array(emity_xs2)

kurt_x_rf  = np.array(kurt_x_rf)
kurt_xp_rf = np.array(kurt_xp_rf)
kurt_x_xs = np.array(kurt_x_xs)
kurt_xp_xs = np.array(kurt_xp_xs)
kurt_x_xs2 = np.array(kurt_x_xs2)
kurt_xp_xs2 = np.array(kurt_xp_xs2)

kurt_t_rf = np.array(kurt_t_rf)
kurt_E_rf = np.array(kurt_E_rf)
kurt_t_xs = np.array(kurt_t_xs)
kurt_E_xs = np.array(kurt_E_xs)
kurt_t_xs2 = np.array(kurt_t_xs2)
kurt_E_xs2 = np.array(kurt_E_xs2)

window_trans = 41
window_long = 11

def moving_average(x, window):
    if len(x) < window:
        return x
    return np.convolve(x, np.ones(window)/window, mode='valid')

#kurt_x_xs_ma   = moving_average(kurt_x_xs, window_trans)
#kurt_xp_xs_ma  = moving_average(kurt_xp_xs, window_trans)
#kurt_t_xs_ma = moving_average(kurt_t_xs, window_long)
#kurt_E_xs_ma = moving_average(kurt_E_xs, window_long)

#kurt_x_xs2_ma  = moving_average(kurt_x_xs2, window_trans)
#kurt_xp_xs2_ma = moving_average(kurt_xp_xs2, window_trans)
#kurt_t_xs2_ma = moving_average(kurt_t_xs2, window_long)
#kurt_E_xs2_ma = moving_average(kurt_E_xs2, window_long)

kurt_x_rf_ma   = moving_average(kurt_x_rf, window_trans)
kurt_xp_rf_ma  = moving_average(kurt_xp_rf, window_trans)
kurt_t_rf_ma = moving_average(kurt_t_rf, window_long)
kurt_E_rf_ma = moving_average(kurt_E_rf, window_long)

emit4d_rf_ma  = moving_average(emit4d_rf, window_trans)
#emit4d_xs_ma  = moving_average(emit4d_xs, window_trans)
#emit4d_xs2_ma = moving_average(emit4d_xs2, window_trans)
emitx_rf_ma = moving_average(emitx_rf, window_trans)
emity_rf_ma = moving_average(emity_rf, window_trans)

emitz_rf_ma  = moving_average(emitz_rf, window_long)
#emitz_xs_ma  = moving_average(emitz_xs, window_long)
#emitz_xs2_ma = moving_average(emitz_xs2, window_long)

# ============================================================
# PLOT COMPARISON
# ============================================================

turns_6d_rf  = rf_turns
turns_6d_xs  = np.arange(1, len(emit6d_xs)+1)
turns_6d_xs2 = np.arange(1, len(emit6d_xs2)+1)

turns_4d_rf  = rf_turns[:len(emit4d_rf_ma)] if len(emit4d_rf_ma) == len(rf_turns) else rf_turns[window_trans//2:window_trans//2+len(emit4d_rf_ma)]
turns_4d_xs  = np.arange(1, len(emit4d_xs)+1)
turns_4d_xs2 = np.arange(1, len(emit4d_xs2)+1)

turns_z_rf  = rf_turns[:len(emitz_rf_ma)] if len(emitz_rf_ma) == len(rf_turns) else rf_turns[window_long//2:window_long//2+len(emitz_rf_ma)]
turns_z_xs  = np.arange(1, len(emitz_xs)+1)
turns_z_xs2 = np.arange(1, len(emitz_xs2)+1)

turns_x_rf = turns_4d_rf
turns_x_xs = np.arange(1,len(emitx_xs)+1)
turns_x_xs2 = np.arange(1,len(emitx_xs2)+1)

# -----------------------------
# 6D emittance
# -----------------------------
plt.figure(figsize=(8,5))

plt.plot(turns_6d_rf, emit6d_rf, 'r', lw=2, label="RF-Track")
plt.plot(turns_6d_xs*25 - 25, emit6d_xs, 'b', lw=2, label="Zampetakis et al.")
plt.plot(turns_6d_xs2*25 - 25, emit6d_xs2, 'k', lw=2, label="Nagaitsev")
plt.legend(fontsize=14, loc='upper left')
plt.grid(True)
plt.xlabel("Number of turns", fontsize=16)
plt.ylabel(r"$\epsilon_{6D}$ ", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("6d_emittance.png", dpi=300, bbox_inches="tight")
plt.savefig("6d_emittance.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(8,5))
plt.plot(turns_4d_rf, emit4d_rf_ma * 500 / mass, 'r', lw=2, label="RF-Track")
plt.plot(turns_4d_xs*25 -25, emit4d_xs * 500 / mass, 'b', lw=2, label="Zampetakis et al.")
plt.plot(turns_4d_xs2*25 - 25, emit4d_xs2 * 500 / mass, 'k', lw=2, label="Nagaitsev")
plt.legend(fontsize=14, loc='upper left')
plt.grid(True)
plt.xlabel("Number of turns", fontsize=16)
plt.ylabel(r"$\epsilon_{4D}$ [mm *mrad]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("4d_transverse_emittance.png", dpi=300, bbox_inches="tight")
plt.savefig("4d_transverse_emittance.pdf", bbox_inches="tight")
plt.close()


# -----------------------------
# Longitudinal emittance
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(turns_z_rf, emitz_rf_ma, 'r', lw=2, label="RF-Track")
plt.plot(turns_z_xs*25 - 25, emitz_xs, 'b', lw=2, label="Zampetakis et al.")
plt.plot(turns_z_xs2*25 - 25, emitz_xs2, 'k', lw=2, label="Nagaitsev")
plt.legend(fontsize=14, loc='upper left')
plt.grid(True)
plt.xlabel("Number of turns", fontsize=16)
plt.ylabel(r"$\epsilon_{z}$ [eV * s]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("longitudinal_emittance.png", dpi=300, bbox_inches="tight")
plt.savefig("longitudinal_emittance.pdf", bbox_inches="tight")
plt.close()

########## TEST ##########
turns_kurt_x_rf  = turns_4d_rf
turns_kurt_x_xs  = np.arange(1, len(kurt_x_xs)+1)
turns_kurt_x_xs2 = np.arange(1, len(kurt_x_xs2)+1)

turns_kurt_z_rf  = turns_z_rf
turns_kurt_z_xs  = np.arange(1, len(kurt_t_xs)+1)
turns_kurt_z_xs2 = np.arange(1, len(kurt_t_xs2)+1)

plt.figure(figsize=(8,5))
plt.plot(turns_kurt_x_rf, kurt_x_rf_ma,'r-', lw=2, label='RF-Track - x')
plt.plot(turns_kurt_x_rf, kurt_xp_rf_ma,'r--', lw=2, label="RF-Track - x'")
plt.plot(turns_kurt_x_xs * 25 - 25, kurt_x_xs,'b-', lw=2, label='Zampetakis et al. (Xsuite) x')
plt.plot(turns_kurt_x_xs * 25 - 25, kurt_xp_xs,'b--', lw=2, label="Zampetakis et al. (Xsuite) x'")
plt.plot(turns_kurt_x_xs2 * 25 -25, kurt_x_xs2,'k-', lw=2, label='Nagaitsev (Xsuite) - x')
plt.plot(turns_kurt_x_xs2 * 25 - 25, kurt_xp_xs2,'k--', lw=2, label="Nagaitsev (Xsuite) - x'")
plt.grid(True)
plt.xlabel("Turn", fontsize=16)
plt.ylabel("Excess kurtosis", fontsize=16)
plt.legend(fontsize=11, ncol=2)
plt.tight_layout()
plt.savefig("kurtosis_transverse.png", dpi=300)
plt.savefig("kurtosis_transverse.pdf")
plt.close()

plt.figure(figsize=(8,5))

plt.plot(turns_kurt_z_rf, kurt_t_rf_ma,'r-', lw=2, label='RF-Track - t')
plt.plot(turns_kurt_z_rf, kurt_E_rf_ma,'r--', lw=2, label='RF-Track - E')
plt.plot(turns_kurt_z_xs * 25 -25, kurt_t_xs,'b-', lw=2, label='Zampetakis et al. (Xsuite) - t')
plt.plot(turns_kurt_z_xs * 25 -25, kurt_E_xs,'b--', lw=2, label='Zampetakis et al. (Xsuite) - E')
plt.plot(turns_kurt_z_xs2 * 25 - 25, kurt_t_xs2,'k-', lw=2, label='Nagaitsev (Xsuite) - t')
plt.plot(turns_kurt_z_xs2 * 25 - 25, kurt_E_xs2,'k--', lw=2, label='Nagaitsev (Xsuite) - E')
plt.grid(True)
plt.xlabel("Turn", fontsize=16)
plt.ylabel("Excess kurtosis", fontsize=16)
plt.legend(fontsize=11, ncol=2)
plt.tight_layout()
plt.savefig("kurtosis_longitudinal.png", dpi=300)
plt.savefig("kurtosis_longitudinal.pdf")
plt.close()

n_total = len(x_rft)

n_trans = np.sum(mask_t_rft)
n_long  = np.sum(mask_z_rft)
n_6d    = np.sum(mask_6_rft)

print("RF-Track:")
print(f"Total particles       : {n_total}")
print(f"Transverse retained   : {n_trans} ({100*n_trans/n_total:.2f}%)")
print(f"Longitudinal retained : {n_long} ({100*n_long/n_total:.2f}%)")
print(f"6D retained           : {n_6d} ({100*n_6d/n_total:.2f}%)")
print(f"Removed particles     : {n_total-n_6d} ({100*(n_total-n_6d)/n_total:.2f}%)")
print("--------------------------------")


plt.figure(figsize=(8,5))
plt.plot(turns_x_rf, emitx_rf_ma*500/mass,
         'r-', lw=2, label=r"RF-Track $\epsilon_{nx}$")
plt.plot(turns_x_rf, emity_rf_ma*500/mass,
         'r--', lw=2, label=r"RF-Track $\epsilon_{ny}$")
plt.plot(turns_x_xs*25-25, emitx_xs*500/mass,
         'b-', lw=2, label=r"Zampetakis et al. (Xsuite) $\epsilon_{nx}$")
plt.plot(turns_x_xs*25-25, emity_xs*500/mass,
         'b--', lw=2, label=r"Zampetakis et al. (Xsuite) $\epsilon_{ny}$")
plt.plot(turns_x_xs2*25-25, emitx_xs2*500/mass,
         'k-', lw=2, label=r"Nagaitsev (Xsuite) $\epsilon_{nx}$")

plt.plot(turns_x_xs2*25-25, emity_xs2*500/mass,
         'k--', lw=2, label=r"Nagaitsev (Xsuite) $\epsilon_{ny}$")
plt.grid(True)
plt.xlabel("Number of turns", fontsize=16)
plt.ylabel(r"$\epsilon_{n}$ [mm mrad]", fontsize=16)
plt.legend(
    fontsize=12,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.18),
    ncol=2,
    frameon=True
)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("transverse_planes_emittance.png",
            dpi=300, bbox_inches="tight")
plt.savefig("transverse_planes_emittance.pdf",
            bbox_inches="tight")
plt.close()