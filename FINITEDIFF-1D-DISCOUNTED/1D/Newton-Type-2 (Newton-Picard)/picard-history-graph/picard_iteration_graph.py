import matplotlib.pyplot as plt

# Data extracted from the provided convergence history
iters = list(range(1, 12))
rel_err_u = [
    3.1623e-02, 4.5168e-03, 6.7421e-04, 6.2634e-04, 3.3921e-04, 
    1.2156e-04, 3.1538e-05, 6.2270e-06, 1.8646e-06, 9.2608e-07, 3.7402e-07
]
rel_err_m = [
    1.5377e-02, 3.2391e-03, 7.1421e-04, 5.3032e-04, 2.6190e-04, 
    9.3222e-05, 2.6658e-05, 7.7174e-06, 3.2580e-06, 1.3906e-06, 5.3026e-07
]

plt.figure(figsize=(10, 6))
plt.semilogy(iters, rel_err_u, 'r-o', label=r"Relative Error of $u$")
plt.semilogy(iters, rel_err_m, 'b-s', label=r"Relative Error of $m$")

# Add a dashed line at 1e-6 and label it epsilon_P at the left
plt.axhline(y=1e-6, color='black', linestyle='--', alpha=0.6, linewidth=1.2)
# Placing the label at the left end of the line
plt.text(iters[0] - 0.5, 1.2e-6, r'$\epsilon_P = 10^{-6}$', color='black', fontsize=12, fontweight='bold')

plt.title('Picard Iteration History')
plt.xlabel('Picard Iterate')
plt.ylabel('Relative Error (Log Scale)')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.xticks(iters)
plt.xlim(0, 12) # Adjusted to allow the label room on the left
plt.legend()

plt.tight_layout()
plt.savefig('discounted1d_newton3_picard_history.png')
plt.show()
