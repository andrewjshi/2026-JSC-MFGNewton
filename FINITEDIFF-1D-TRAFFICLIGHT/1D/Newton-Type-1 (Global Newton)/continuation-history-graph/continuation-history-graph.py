import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Hardcoded values extracted from newton-type1-history_trafficlight_T1.0_nu0.001.txt
# Format: (viscosity_value, convergence_status)
attempts = [
    (0.1000, True),  # Attempting nu = 0.1000
    (0.0500, True),  # Attempting nu = 0.0500
    (0.0250, True),  # Attempting nu = 0.0250
    (0.0125, True),  # Attempting nu = 0.0125
    (0.0063, False), # Attempting nu = 0.0063 (Failed)
    (0.0094, True),  # Attempting nu = 0.0094
    (0.0047, False), # Attempting nu = 0.0047 (Failed)
    (0.0070, True),  # Attempting nu = 0.0070
    (0.0035, False), # Attempting nu = 0.0035 (Failed)
    (0.0053, True),  # Attempting nu = 0.0053
    (0.0026, False), # Attempting nu = 0.0026 (Failed)
    (0.0040, False), # Attempting nu = 0.0040 (Failed)
    (0.0046, True),  # Attempting nu = 0.0046
    (0.0023, False), # Attempting nu = 0.0023 (Failed)
    (0.0035, False), # Attempting nu = 0.0035 (Failed)
    (0.0040, True),  # Attempting nu = 0.0040
    (0.0020, False), # Attempting nu = 0.0020 (Failed)
    (0.0030, False), # Attempting nu = 0.0030 (Failed)
    (0.0035, False), # Attempting nu = 0.0035 (Failed)
    (0.0038, False), # Attempting nu = 0.0038 (Failed)
    (0.0039, True),  # Attempting nu = 0.0039
    (0.0020, False), # Attempting nu = 0.0020 (Failed)
    (0.0029, False), # Attempting nu = 0.0029 (Failed)
    (0.0034, False), # Attempting nu = 0.0034 (Failed)
    (0.0037, False), # Attempting nu = 0.0037 (Failed)
    (0.0038, True),  # Attempting nu = 0.0038
    (0.0019, False), # Attempting nu = 0.0019 (Failed)
    (0.0028, False), # Attempting nu = 0.0028 (Failed)
    (0.0033, False), # Attempting nu = 0.0033 (Failed)
    (0.0036, False), # Attempting nu = 0.0036 (Failed)
    (0.0037, True),  # Attempting nu = 0.0037
    (0.0018, False), # Attempting nu = 0.0018 (Failed)
    (0.0028, False), # Attempting nu = 0.0028 (Failed)
    (0.0032, False), # Attempting nu = 0.0032 (Failed)
    (0.0034, True),  # Attempting nu = 0.0034
    (0.0017, False), # Attempting nu = 0.0017 (Failed)
    (0.0026, False), # Attempting nu = 0.0026 (Failed)
    (0.0030, False), # Attempting nu = 0.0030 (Failed)
    (0.0032, False), # Attempting nu = 0.0032 (Failed)
    (0.0033, True),  # Attempting nu = 0.0033
    (0.0017, False), # Attempting nu = 0.0017 (Failed)
    (0.0025, False), # Attempting nu = 0.0025 (Failed)
    (0.0029, False), # Attempting nu = 0.0029 (Failed)
    (0.0031, False), # Attempting nu = 0.0031 (Failed)
    (0.0032, True),  # Attempting nu = 0.0032
    (0.0016, False), # Attempting nu = 0.0016 (Failed)
    (0.0024, False), # Attempting nu = 0.0024 (Failed)
    (0.0028, False), # Attempting nu = 0.0028 (Failed)
    (0.0030, False), # Attempting nu = 0.0030 (Failed)
    (0.0031, False), # Attempting nu = 0.0031 (Failed)
]

target_nu = 0.001

# Generate the plot
x = range(1, len(attempts) + 1)
y = [a[0] for a in attempts]
colors = ['black' if a[1] else 'red' for a in attempts]

plt.figure(figsize=(10, 6))

# Plot connecting line
plt.plot(x, y, color='gray', alpha=0.5, linestyle='-', zorder=1)

# Plot the dots
plt.scatter(x, y, c=colors, s=60, zorder=2) # Slightly smaller dots since there are 51 points

# Target line
plt.axhline(y=target_nu, color='black', linestyle=':', linewidth=2, zorder=1)

# Formatting - Adjusted for the new nu range (0 to 0.20)
plt.ylim(-0.01, 0.12)
plt.yticks([0.0, 0.05, 0.10])

plt.xlabel('Continuation Attempt')
plt.ylabel(r'$\nu$ (Viscosity)')
plt.title('Viscosity Continuation Progress (Newton Type 1 - Traffic Light)')
plt.grid(True, linestyle='--', alpha=0.6)

# Reduce x-ticks crowding by showing every 5th tick
plt.xticks(range(1, len(attempts) + 1, 5))

# Custom legend
custom_lines = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10),
    Line2D([0], [0], color='black', linestyle=':', linewidth=2)
]
plt.legend(custom_lines, ['Converged', 'Failed', r'Target $\nu = 0.001$'], loc='upper right')

plt.tight_layout()
plt.savefig('trafficlight_newton1_history.png')
plt.show()

