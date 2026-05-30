import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Hardcoded values extracted from newton-type1-history_discounted_T10.0_nu0.5.txt
# Format: (viscosity_value, convergence_status)
attempts = [
    (1.0, True),     # First Solve: nu = 1.0000 (Converged)
    (0.5, False),    # Attempting nu = 0.5000 (Failed)
    (0.75, True),    # Attempting nu = 0.7500 (Converged)
    (0.5, False),    # Attempting nu = 0.5000 (Failed)
    (0.625, True),   # Attempting nu = 0.6250 (Converged)
    (0.5, False),    # Attempting nu = 0.5000 (Failed)
    (0.5625, True),  # Attempting nu = 0.5625 (Converged)
    (0.5, True)      # Attempting nu = 0.5000 (Converged)
]

target_nu = 0.5

# Generate the plot
x = range(1, len(attempts) + 1)
y = [a[0] for a in attempts]
colors = ['black' if a[1] else 'red' for a in attempts]

plt.figure(figsize=(10, 6))

# Plot connecting line
plt.plot(x, y, color='gray', alpha=0.5, linestyle='-', zorder=1)

# Plot the dots
plt.scatter(x, y, c=colors, s=100, zorder=2)

# Target line
plt.axhline(y=target_nu, color='black', linestyle=':', linewidth=2, zorder=1)

# Formatting
plt.ylim(0.0, 1.1)
plt.yticks([0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0])

plt.xlabel('Continuation Attempt')
plt.ylabel(r'$\nu$ (Viscosity)')
plt.title('Viscosity Continuation Progress (Newton Type 1)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(x)

# Custom legend
custom_lines = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10),
    Line2D([0], [0], color='black', linestyle=':', linewidth=2)
]
plt.legend(custom_lines, ['Converged', 'Failed', r'Target $\nu = 0.5$'], loc='upper right')

plt.tight_layout()
plt.savefig('discounted_newton1_newton_history.png')
plt.show()