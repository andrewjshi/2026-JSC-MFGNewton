import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Values extracted from newton-type1-history_doublewell_T1.0_nu0.01.txt
# Format: (viscosity_value, convergence_status)
attempts = [
    (0.5, True), (0.375, False), (0.4375, True), (0.3281, False), (0.3828, True), 
    (0.2871, False), (0.335, True), (0.2512, False), (0.2931, True), (0.2198, False), 
    (0.2565, True), (0.1923, False), (0.2244, True), (0.1683, True), (0.1262, False), 
    (0.1473, True), (0.1104, False), (0.1289, True), (0.0966, False), (0.1127, True), 
    (0.0846, False), (0.0987, True), (0.074, False), (0.0863, True), (0.0647, False), 
    (0.0755, True), (0.0566, False), (0.0661, True), (0.0496, False), (0.0578, True), 
    (0.0434, False), (0.0506, True), (0.038, False), (0.0443, True), (0.0332, False), 
    (0.0387, True), (0.0291, False), (0.0339, True), (0.0254, False), (0.0297, True), 
    (0.0222, False), (0.026, True), (0.0195, False), (0.0227, True), (0.017, True), 
    (0.0128, True), (0.01, True)
]

target_nu = 0.01

# Generate the plot
x = range(1, len(attempts) + 1)
y = [a[0] for a in attempts]
colors = ['black' if a[1] else 'red' for a in attempts]

plt.figure(figsize=(12, 6))

# Plot connecting line
plt.plot(x, y, color='gray', alpha=0.5, linestyle='-', zorder=1)

# Plot the dots
plt.scatter(x, y, c=colors, s=60, zorder=2)

# Target line
plt.axhline(y=target_nu, color='black', linestyle=':', linewidth=2, zorder=1)

# Formatting
plt.xlabel('Continuation Attempt')
plt.ylabel(r'$\nu$ (Viscosity)')
plt.title('Viscosity Continuation Progress (Newton Type 1 - Double Well)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(x[::2]) # Subsample x-ticks for better readability

# Custom legend
custom_lines = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10),
    Line2D([0], [0], color='black', linestyle=':', linewidth=2)
]
plt.legend(custom_lines, ['Converged', 'Failed', r'Target $\nu = 0.01$'], loc='upper right')

plt.tight_layout()
plt.savefig('doublewell_newton1_continuation_history.png')
plt.show()

