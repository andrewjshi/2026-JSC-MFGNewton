import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import time
import sys

# ==========================================
# SECTION 1: PARAMETERS AND GRID DISCRETIZATION
# ==========================================

# --- Physics Parameters ---
rho = 0.5                # Discount factor 
mu = 0                   # Analytical center

# --- Spatial and Temporal Grid ---
width = 2.0              
xmin, xmax = mu - width, mu + width # Spatial domain boundaries
T = 10.0                 # Total time horizon
Nx = 200                 # Number of spatial intervals
Nt = 200                 # Number of temporal intervals

# --- Newton Iteration Parameters ---
NEWTON_MAX_ITER = 20     # Maximum Newton iterations for the backward HJB solves
NEWTON_TOL = 1e-9        # Convergence tolerance for the backward Newton solver

# --- Continuation Parameters ---
NU_START = 1.0           # Starting viscosity
NU_TARGET = 0.5          # Target viscosity
BETA = 0.5               # Reduction factor for viscosity step
GAMMA = 0.5              # Backtracking factor
CONT_TOL = 1e-4          # Minimum allowable step size

# --- Derived Quantities ---
L = xmax - xmin          # Length of the spatial domain
Dx = L / Nx              # Spatial step size
dt = T / Nt              # Temporal step size            
x = np.linspace(xmin, xmax, Nx + 1) # Computational grid points (Dirichlet)
num_x = x.size

# --- Saved File Names ---
history_filename = f"newton-type1-history_discounted1d_T{T}_nu{NU_TARGET}.txt"
plot_snapshot_filename = f"type1-discounted1d_plot_T{T}_nu{NU_TARGET}.png"

# ==========================================
# SECTION 2: PHYSICS DEFINITIONS (INCLUDING HAMILTONIAN) & EXACT SOLUTION
# ==========================================

def get_exact_params(nu_val):
    eta = (1.0 + rho * nu_val) / (2.0 * nu_val)
    s_sq = nu_val / (2.0 * eta)
    m_peak = 1.0 / np.sqrt(2.0 * np.pi * s_sq)
    omega = (2.0 * nu_val * eta + np.log(m_peak)) / rho
    return eta, s_sq, m_peak, omega

def get_exact_m(x_vals, nu_val):
    _, s_sq, m_peak, _ = get_exact_params(nu_val)
    return m_peak * np.exp(-(x_vals - mu)**2 / (2.0 * s_sq))

def get_exact_u(x_vals, nu_val):
    eta, _, _, omega = get_exact_params(nu_val)
    return eta * (x_vals - mu)**2 + omega

m_star = get_exact_m(x, NU_TARGET)
u_star = get_exact_u(x, NU_TARGET)

def compute_f(m_dist):
    m_safe = np.maximum(m_dist, 1e-12)
    return -np.log(m_safe)

def discrete_Hamiltonian(dp, dm, f_val):
    H_p = 0.5 * np.minimum(dp, 0)**2 + 0.5 * np.maximum(dm, 0)**2
    return f_val - H_p 

def dH_dp(p):
    return p

# ==========================================
# SECTION 3: DISCRETE OPERATORS
# ==========================================

I = sp.eye(num_x)

# Forward Difference Matrix (D+) - Dirichlet
data_p = np.array([-np.ones(num_x), np.ones(num_x)])
D_plus = sp.spdiags(data_p, [0, 1], num_x, num_x).tolil()
D_plus[num_x-1, :] = 0 
D_plus = D_plus.tocsr() / Dx

# Backward Difference Matrix (D-) - Dirichlet
data_m = np.array([-np.ones(num_x), np.ones(num_x)])
D_minus = sp.spdiags(data_m, [-1, 0], num_x, num_x).tolil()
D_minus[0, :] = 0 
D_minus = D_minus.tocsr() / Dx

# Laplacian Matrix
Laplacian = D_plus @ D_minus

# ==========================================
# SECTION 4: GLOBAL NEWTON METHOD (NEWTON TYPE 1) & VISCOSITY CONTINUATION
# ==========================================

def newton_type1_step(nu_val, W_guess, log_file=None):
    m_st, u_st = get_exact_m(x, nu_val), get_exact_u(x, nu_val)
    m_bnd, u_bnd = m_st[0], u_st[0]
    
    if W_guess is None:
        U = np.full((Nt + 1, num_x), u_bnd)
        M = np.full((Nt + 1, num_x), m_bnd)
    else:
        U = W_guess[:(Nt+1)*num_x].reshape((Nt+1, num_x))
        M = W_guess[(Nt+1)*num_x:].reshape((Nt+1, num_x))

    I_neg = -sp.eye(num_x).tolil()
    I_neg[0, 0] = 0; I_neg[-1, -1] = 0
    I_neg = I_neg.tocsr()

    dv_du = -0.5 * (D_plus + D_minus)

    for it in range(NEWTON_MAX_ITER):
        J_UU = [[None for _ in range(Nt+1)] for _ in range(Nt+1)]
        J_UM = [[None for _ in range(Nt+1)] for _ in range(Nt+1)]
        J_MU = [[None for _ in range(Nt+1)] for _ in range(Nt+1)]
        J_MM = [[None for _ in range(Nt+1)] for _ in range(Nt+1)]

        F_U, F_M = np.zeros((Nt+1, num_x)), np.zeros((Nt+1, num_x))
        A_diff_U = (1.0 + rho * dt) * I - dt * nu_val * Laplacian

        J_UU[Nt][Nt], F_U[Nt] = sp.eye(num_x).tocsr(), U[Nt] - u_st
        J_MM[0][0], F_M[0] = sp.eye(num_x).tocsr(), M[0] - m_st
        J_UM[Nt][Nt] = J_MU[0][0] = J_MU[Nt][Nt] = sp.csr_matrix((num_x, num_x))

        for n in range(Nt):
            u_curr, u_next, m_curr, m_next = U[n], U[n+1], M[n], M[n+1]
            dp, dm = D_plus @ u_curr, D_minus @ u_curr
            p_min, p_max, m_safe = np.minimum(dp, 0), np.maximum(dm, 0), np.maximum(m_curr, 1e-12)
            
            F_u = A_diff_U @ u_curr + dt * discrete_Hamiltonian(dp, dm, compute_f(m_safe)) - u_next
            F_u[0], F_u[-1] = u_curr[0] - u_bnd, u_curr[-1] - u_bnd
            F_U[n] = F_u

            J_u = A_diff_U - dt * (sp.diags(dH_dp(p_min)) @ D_plus + sp.diags(dH_dp(p_max)) @ D_minus)
            J_u = J_u.tolil(); J_u[0,:], J_u[0,0], J_u[-1,:], J_u[-1,-1] = 0, 1, 0, 1
            J_UU[n][n], J_UU[n][n+1] = J_u.tocsr(), I_neg
            J_um = -dt * sp.diags(1.0 / m_safe).tolil(); J_um[0,:], J_um[-1,:] = 0, 0
            J_UM[n][n] = J_um.tocsr()

            v = -dH_dp((dp + dm) / 2.0) 
            v_p, v_m = np.maximum(v, 0), np.minimum(v, 0)
            Adv_Op = D_minus @ sp.diags(v_p) + D_plus @ sp.diags(v_m)
            A_M = I - dt * nu_val * Laplacian + dt * Adv_Op

            F_m = A_M @ m_next - m_curr
            F_m[0], F_m[-1] = m_next[0] - m_bnd, m_next[-1] - m_bnd
            F_M[n+1] = F_m

            A_bnd_M = A_M.tolil(); A_bnd_M[0,:], A_bnd_M[0,0], A_bnd_M[-1,:], A_bnd_M[-1,-1] = 0, 1, 0, 1
            J_MM[n+1][n+1], J_MM[n+1][n] = A_bnd_M.tocsr(), I_neg
            
            Hv = (v > 0).astype(float)
            dAdv_dv = D_minus @ sp.diags(m_next * Hv) + D_plus @ sp.diags(m_next * (1.0 - Hv))
            J_mu = dt * (dAdv_dv.tocsr() @ dv_du).tolil(); J_mu[0,:], J_mu[-1,:] = 0, 0
            J_MU[n+1][n] = J_mu.tocsr()

        J = sp.bmat([[sp.bmat(J_UU), sp.bmat(J_UM)], [sp.bmat(J_MU), sp.bmat(J_MM)]]).tocsr()
        F = np.concatenate([F_U.flatten(), F_M.flatten()])
        res = np.linalg.norm(F, np.inf)
        
        if np.isnan(res) or np.isinf(res): return np.concatenate([U.flatten(), M.flatten()]), False, res
        if res < NEWTON_TOL: 
            iter_log = f"      Iter {it:2d} | res = {res:.4e} (Converged)"
            print(iter_log)
            if log_file: log_file.write(iter_log + "\n")
            return np.concatenate([U.flatten(), M.flatten()]), True, res

        # --- TIMER FOR LINEAR SOLVE ---
        t_solve_start = time.time()
        try:
            dW = spla.spsolve(J, F)
            t_solve_elapsed = time.time() - t_solve_start
            
            iter_log = f"      Iter {it:2d} | res = {res:.4e} | Solve Time: {t_solve_elapsed:.4f}s"
            print(iter_log)
            if log_file: log_file.write(iter_log + "\n")

            U = (U.flatten() - dW[:(Nt+1)*num_x]).reshape((Nt+1, num_x))
            M = (M.flatten() - dW[(Nt+1)*num_x:]).reshape((Nt+1, num_x))
        except: 
            return np.concatenate([U.flatten(), M.flatten()]), False, res

    return np.concatenate([U.flatten(), M.flatten()]), False, res

# --- Viscosity Continuation Main Loop ---
start_time_all = time.time()
with open(history_filename, "w") as f:
    header = (
        f"{'='*118}\n"
        f"   DISCOUNTED MFG EXAMPLE (NEWTON TYPE 1: GLOBAL NEWTON)\n"
        f"{'='*118}\n"
        f"Parameters:\n"
        f"  xmin = {xmin}, xmax = {xmax}\n"
        f"  T    = {T}, Nx = {Nx}, Nt = {Nt}\n"
        f"  nu   = {NU_TARGET}, rho = {rho}, mu = {mu}\n"
        f"Solver Parameters:\n"
        f"  NEWTON_MAX_ITER = {NEWTON_MAX_ITER}, NEWTON_TOL = {NEWTON_TOL}\n"
        f"  Continuation: NU_START={NU_START}, BETA={BETA}, GAMMA={GAMMA}\n"
        f"Grid Info:\n"
        f"  dt = {dt:.6f}, dx = {Dx:.6f}\n"
        f"{'-'*118}\n"
    )
    print(header, end=""); f.write(header)
    
    nu_curr = NU_START
    print(f"First Solve: nu = {nu_curr:.4f}"); f.write(f"First Solve: nu = {nu_curr:.4f}\n")
    
    t_nu_start = time.time()
    W, converged, res = newton_type1_step(nu_curr, None, f)
    t_nu_elapsed = time.time() - t_nu_start
    
    if not converged: sys.exit()
    conv_msg = f"    -> Converged! (res={res:.4e}) Total for nu={nu_curr:.4f}: {t_nu_elapsed:.4f}s\n"
    print(conv_msg); f.write(conv_msg + "\n")

    while nu_curr > NU_TARGET:
        nu_next = max(NU_TARGET, BETA * nu_curr)
        success = False
        while not success:
            print(f"  Attempting nu = {nu_next:.4f}..."); f.write(f"  Attempting nu = {nu_next:.4f}...\n")
            
            t_nu_start = time.time()
            W_trial, conv, r = newton_type1_step(nu_next, W, f)
            t_nu_elapsed = time.time() - t_nu_start
            
            if conv:
                W, nu_curr, success = W_trial, nu_next, True
                conv_msg = f"    -> Converged! (res={r:.4e}) Total for nu={nu_next:.4f}: {t_nu_elapsed:.4f}s\n"
                print(conv_msg); f.write(conv_msg + "\n")
            else:
                fail_msg = f"    -> Failed (res={r:.4e}) in {t_nu_elapsed:.4f}s. Backtracking...\n"
                print(fail_msg); f.write(fail_msg)
                nu_next = nu_curr - GAMMA * (nu_curr - nu_next)
                if abs(nu_curr - nu_next) < CONT_TOL: sys.exit()

    total_time = time.time() - start_time_all
    time_msg = f"{'-'*118}\nTotal Execution Time: {total_time:.4f} seconds.\n"
    print(time_msg); f.write(time_msg)

# ==========================================
# SECTION 5: VISUALIZATION
# ==========================================

U_flow = W[:(Nt+1)*num_x].reshape((Nt+1, num_x))
M_flow = W[(Nt+1)*num_x:].reshape((Nt+1, num_x))

idx_mid = Nt // 2
fig2 = plt.figure(figsize=(14, 5))
plt.subplot(1, 2, 1)
plt.plot(x, u_star, 'k--', linewidth=2, label='Stationary analytical $u^*$')
plt.plot(x, U_flow[idx_mid], 'r-', linewidth=2, label=f'Numerical $u(T/2)$')
plt.title(f'Value Function Comparison (T={T})')
plt.xlabel('x'); plt.legend(); plt.grid(True)
plt.subplot(1, 2, 2)
plt.plot(x, m_star, 'k--', linewidth=2, label='Stationary analytical $m^*$')
plt.plot(x, M_flow[idx_mid], 'b-', linewidth=2, label=f'Numerical $m(T/2)$')
plt.title(f'Density Comparison (T={T})')
plt.xlabel('x'); plt.legend(); plt.grid(True)
plt.tight_layout(); plt.savefig(plot_snapshot_filename); plt.close(fig2)

# Save the 1D arrays to data files for TikZ/LaTeX plotting
np.savetxt('type1_discounted_value_function_comparison.dat', np.column_stack((x, u_star, U_flow[idx_mid])), header='x u_star u_num', comments='', fmt='%.6e')
np.savetxt('type1_discounted_density_comparison.dat', np.column_stack((x, m_star, M_flow[idx_mid])), header='x m_star m_num', comments='', fmt='%.6e')    

    