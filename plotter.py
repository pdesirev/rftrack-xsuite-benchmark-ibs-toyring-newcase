import numpy as np
import matplotlib.pyplot as plt

# Load data
rft = np.loadtxt('Results/RFT_25pC_nemit_x_IBS_1000.dat')  # already in mm·mrad
xsuite = np.loadtxt('Results/XSuite_25pC_nemitt_x_1000.dat')    # in m·rad
xsuite_analytical = np.loadtxt('Results/XSuite_analytical_25pC_nemitt_x_1000.dat')

rft_y = np.loadtxt('Results/RFT_25pC_nemit_y_IBS_1000.dat')
xsuite_y = np.loadtxt('Results/XSuite_25pC_nemitt_y_1000.dat')
xsuite_analytical_y = np.loadtxt('Results/XSuite_analytical_25pC_nemitt_y_1000.dat')

# Convert XSuite from m·rad → mm·mrad
xsuite *= 1e6
xsuite_analytical *= 1e6
xsuite_y *= 1e6
xsuite_analytical_y *= 1e6

# x-axis
x = np.arange(len(rft))

# Plot
plt.figure(figsize=(10, 6))

# Horizontal (x)
plt.plot(x, rft, 'r', linewidth=2.0, label=r'RF-Track $\epsilon_x$')
plt.plot(x, xsuite, 'b', linewidth=2.0, label=r'Zampetakis et al. $\epsilon_x$')
plt.plot(x, xsuite_analytical, 'k', linewidth=2.0, label=r'Nagaitsev $\epsilon_x$')

# Vertical (y)
plt.plot(x, rft_y, 'b--', linewidth=2.0, label=r'RF-Track $\epsilon_y$')
plt.plot(x, xsuite_y, 'r--', linewidth=2.0, label=r'Zampetakis et al. $\epsilon_y$')
plt.plot(x, xsuite_analytical_y, 'k--', linewidth=2.0, label=r'Nagaitsev $\epsilon_y$')

# Labels and style
plt.xlabel('Number of turns', fontsize=18)
plt.ylabel('Normalized emittance [mm·mrad]', fontsize=18)

plt.title('25 pC, IBS', fontsize=18)

plt.legend(fontsize=14, loc='upper left')
plt.grid(True)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Save
plt.savefig('normalized_emittance.png', dpi=300, bbox_inches='tight')
plt.savefig('normalized_emittance.pdf', bbox_inches='tight')

plt.show()
