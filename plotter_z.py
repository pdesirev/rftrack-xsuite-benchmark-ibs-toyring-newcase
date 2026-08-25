import numpy as np
import matplotlib.pyplot as plt
import os

mass = 0.51099895069
Q = -1
pop = 25  # RF_Track.pC not needed numerically here

sigma_z = np.zeros(1000)
emitt_z = np.zeros(1000)

cut = 3

for i in range(1, 1001):
    filename = f"Results/ps/phase_space_{i}.dat"
    
    if not os.path.exists(filename):
        print(f"File missing: {filename}")
        continue

    M = np.loadtxt(filename)
    M = M.reshape(-1, 6)

    P = M[:, 5]  # MATLAB col 6 → Python index 5
    E = np.sqrt(P**2 + mass**2) * 1e6
    mE = np.mean(E)

    gamma = mE / mass
    beta = np.sqrt(1.0 - 1.0 / (gamma**2))

    t = M[:, 4] * 1e-03 * beta  # MATLAB col 5 → index 4
    mP = np.mean(P)

    delta = (P - mP) / mP

    mt = np.mean(t)
    st = np.std(t)
    md = np.mean(delta)
    sd = np.std(delta)

    mask = (np.abs(t - mt) < cut * st) & (np.abs(delta - md) < cut * sd)

    Z = np.column_stack((t[mask], delta[mask]))

    if Z.shape[0] > 1:
        detZ = np.linalg.det(np.cov(Z, rowvar=False))
        emitt_z[i - 1] = np.sqrt(abs(detZ))
    else:
        emitt_z[i - 1] = np.nan  # safety

# x-axis
x = np.arange(1, 1001)

# Load comparison data
XSuite_gemittz_a = np.loadtxt('Results/XSuite_analytical_25pC_gemitt_z_1000.dat')
XSuite_gemittz_m = np.loadtxt('Results/XSuite_25pC_gemitt_z_1000.dat')

# Plot
plt.figure(figsize=(10, 6))

plt.plot(x, emitt_z, 'r', linewidth=2.0, label='RF-Track')
plt.plot(x, XSuite_gemittz_m, 'b', linewidth=2.0, label='Zampetakis et al.')
plt.plot(x, XSuite_gemittz_a, 'k', linewidth=2.0, label='Nagaistev')

plt.legend(fontsize=14, loc='upper left')
plt.grid(True)

plt.xlabel('Number of turns', fontsize=16)
plt.ylabel('Longitudinal emittance [m]', fontsize=16)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Save figure
plt.savefig('longitudinal_emittance.png', dpi=300, bbox_inches='tight')
plt.savefig('longitudinal_emittance.pdf', bbox_inches='tight')

plt.show()
