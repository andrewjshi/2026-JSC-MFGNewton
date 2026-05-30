import matplotlib.pyplot as plt

# Data extracted from newton-convergence-history_doublewell_T1.0_nu0.01.txt
iters = list(range(1, 15))

# Relative Error of U (Column 3)
rel_err_u = [
    6.3679e+09, 3.1298e-03, 4.8646e-04, 9.7534e-05, 7.9960e-05, 
    7.6616e-05, 6.8954e-05, 6.0895e-05, 5.0712e-05, 3.4362e-05, 
    1.2789e-05, 6.0066e-06, 2.2690e-06, 5.1292e-07
]

# Relative Error of M (Column 5)
rel_err_m = [
    2.0282e-02, 3.7597e-03, 1.0858e-03, 1.6128e-04, 2.9340e-05, 
    8.4523e-06, 1.8322e-06, 4.4196e-07, 9.1577e-08, 2.3410e-08, 
    5.3675e-09, 1.2856e-09, 3.6248e-10, 8.4300e-11
]

plt.figure(figsize=(10, 6))

# Plotting u starting from the second iterate (index 1) to exclude the outlier
plt.semilogy(iters[1:], rel_err_u[1:], 'r-o', label=r"Relative Error of $u$")

# Plotting m starting from the first iterate
plt.semilogy(iters, rel_err_m, 'b-s', label=r"Relative Error of $m$")

# Add a dashed line at 1e-6 and label it epsilon_P
plt.axhline(y=1e-6, color='black', linestyle='--', alpha=0.6, linewidth=1.2)
plt.text(iters[0] + 0.1, 1.2e-6, r'$\epsilon_P = 10^{-6}$', color='black', fontsize=12, fontweight='bold')

plt.title('Picard Iteration History (Double Well Example)')
plt.xlabel('Picard Iterate')
plt.ylabel('Relative Error (Log Scale)')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.xticks(iters)
plt.xlim(0.5, 14.5)
plt.legend()

plt.tight_layout()
plt.savefig('doublewell_newton3_picard_history.png')
plt.show()

