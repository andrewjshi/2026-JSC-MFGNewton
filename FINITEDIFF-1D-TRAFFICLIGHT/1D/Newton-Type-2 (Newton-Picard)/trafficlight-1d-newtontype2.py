"""
=============================================================================
PROBLEM SUMMARY: 1D Traffic Light Mean Field Game (Hard Walls) - High Difficulty
=============================================================================

1. SCENARIO
   - A crowd of agents starts at x=0.1 and wants to reach x=1.0.
   - A "Traffic Light" at x=0.5 turns RED during t=[0.3, 0.7].
   - Agents must optimize their speed to reach the target while minimizing
     effort and avoiding the high cost of the red light.

2. PDE SYSTEM (Mean Field Game Equations)
   The problem solves the coupled system for the Value Function u(x,t) 
   and the Density m(x,t):

   (A) Backward Hamilton-Jacobi-Bellman (HJB) Equation:
       Describes the optimal control strategy.
       
       -∂_t u - ν ∂_xx u + ½|∂_x u|² = V(x,t) + m^β
       
       Where:
       - ν (nu): Viscosity/Noise parameter.
       - ½|∂_x u|²: The Hamiltonian H(p) representing kinetic energy cost.
       - V(x,t): The external potential (The Traffic Light).
       - m^β: The congestion cost (agents dislike crowds).

   (B) Forward Fokker-Planck (FP) Equation:
       Describes the evolution of the crowd density under optimal control.
       
       ∂_t m - ν ∂_xx m + ∂_x (m · v) = 0
       
       Where:
       - v(x,t) = -∂_x u  (The optimal velocity field).

3. CONDITIONS
   
   (A) Boundary Conditions (Neumann / Hard Walls):
       The domain x ∈ [0, 1] is closed. Agents cannot leave or loop.
       - ∂_x u(0,t) = ∂_x u(1,t) = 0  (Zero velocity at walls)
       - ∂_x m(0,t) = ∂_x m(1,t) = 0  (Zero flux at walls)

   (B) Initial Condition (t=0):
       - m(x,0) = m0(x) 
       - A Gaussian distribution centered at x=0.1 (Traffic Jam start).

   (C) Terminal Condition (t=T):
       - u(x,T) = G(x) = 5.0 * (x - 1.0)²
       - A quadratic penalty that forces agents to be near x=1.0 at the end.

4. THE TRAFFIC LIGHT (V(x,t))
   - Modeled as a time-dependent Gaussian barrier at x=0.5.
   - V(x,t) = 0 when Green.
   - V(x,t) = High Cost when Red (t ∈ [0.3, 0.7]).
   - Activation is instantaneous (step function).

=============================================================================
"""

import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

# ==========================================
# SECTION 1: PARAMETERS AND GRID DISCRETIZATION
# ==========================================

# --- Physics Parameters ---
nu = 0.02                # Viscosity coefficient
beta = 2                 # Congestion exponent

# --- Traffic Light (Instantaneous) ---
LIGHT_POS = 0.5          # Position of the light
LIGHT_START = 0.3        # Time light turns RED
LIGHT_END = 0.7          # Time light turns GREEN
LIGHT_PENALTY = 150.0    # Cost of running red light (Extreme barrier)
LIGHT_WIDTH = 0.05       # Width of the stop line

# --- Spatial and Temporal Grid ---
xmin, xmax = 0.0, 1.0    # Spatial domain boundaries
T  = 1.0                 # Total time horizon
Nx = 200                 # Number of spatial intervals
Nt = 100                 # Number of temporal intervals

# --- Picard and Newton Iteration Parameters ---
PICARD_MAX_ITER = 200    # Maximum number of outer HJB-FP coupling iterations 
PICARD_TOL = 1e-6        # Convergence tolerance for the outer Picard loop
DAMPING = 0.1            # Damping parameter (0 < DAMPING <= 1). 1 corresponds to no damping. (Heavy damping to prevent bang-bang oscillation)

NEWTON_MAX_ITER = 20     # Maximum Newton iterations for the backward HJB solves
NEWTON_TOL = 1e-9        # Convergence tolerance for the backward Newton solver

# --- Derived Quantities ---
L  = xmax - xmin         # Length of the spatial domain
Dx = L / Nx              # Spatial step size
dt = T / Nt              # Temporal step size      
x  = np.linspace(xmin, xmax, Nx) # Computational grid points (Neumann)
num_x = x.size
norm_const = np.sqrt(Dx * dt)

# --- Directories & Saved File Names ---
IMG_DIR = "img"
os.makedirs(IMG_DIR, exist_ok=True)

history_filename = f"newton-type2-history_trafficlight1d_T{T}_nu{nu}.txt"
plot_snapshot_filename = f"trafficlight1d-plot_T{T}_nu{nu}.png"
video_filename = f"traffic_light_hardwall_T{T}_nu{nu}.mp4"

# ==========================================
# SECTION 2: PHYSICS DEFINITIONS (INCLUDING HAMILTONIAN)
# ==========================================

# --- 1. The Instantaneous "Traffic Light" Potential ---
V_field = np.zeros((Nt + 1, num_x))
barrier_x = np.exp(-((x - LIGHT_POS)**2) / (2 * (LIGHT_WIDTH/2)**2))

for n in range(Nt + 1):
    t = n * dt
    # Instantaneous Step Activation
    intensity = 1.0 if (LIGHT_START <= t <= LIGHT_END) else 0.0
    
    # The Red Light Barrier (Base pull -x removed)
    V_field[n, :] = (LIGHT_PENALTY * intensity * barrier_x)

# --- 2. Initial Config (Traffic Jam on Left) ---
m0 = np.exp(-((x - 0.1)**2) / 0.01)
m0 = m0 / (np.sum(m0) * Dx)

# --- 3. Terminal Cost (Desire to reach Right) ---
# High cost at Left, Zero cost at Right
uT = 5.0 * (x - 1.0)**2 

def discrete_Hamiltonian(p):
    return 0.5 * p**2

def dH_dp(p):
    return p

# ==========================================
# SECTION 3: DISCRETE OPERATORS (NEUMANN / HARD WALLS)
# ==========================================

I = sp.eye(num_x)

# -- LAPLACIAN (D2) -- Nuemann
D2 = sp.spdiags([np.ones(num_x), -2*np.ones(num_x), np.ones(num_x)], [-1, 0, 1], num_x, num_x, format='csr') / (Dx**2)

# Left Wall (x=0): Reflection u_{-1} = u_{1} -> (2u_{1} - 2u_{0})
D2[0, 0] = -2 / (Dx**2)
D2[0, 1] =  2 / (Dx**2)
D2[0, -1] = 0 # Remove periodic wrap

# Right Wall (x=1): Reflection u_{N} = u_{N-2} -> (2u_{N-2} - 2u_{N-1})
D2[-1, -1] = -2 / (Dx**2)
D2[-1, -2] =  2 / (Dx**2)
D2[-1, 0] = 0 # Remove periodic wrap

# -- GRADIENT (Grad) --
# Interior: Central Difference
Grad = sp.spdiags([-0.5*np.ones(num_x), 0.5*np.ones(num_x)], [-1, 1], num_x, num_x, format='csr') / Dx

# Boundaries: Zero Flux (Velocity = 0 at walls)
Grad[0, :] = 0
Grad[-1, :] = 0

# ==========================================
# SECTION 4: HJB AND FP SOLVERS
# ==========================================

def solve_hjb_backward(M_flow):
    U = np.zeros((Nt + 1, Nx))
    U[Nt] = uT 
    
    A_diff = sp.eye(Nx) - dt * nu * D2

    total_newton = 0
    for n in range(Nt - 1, -1, -1):
        u_next = U[n+1]
        u_curr = u_next.copy()
        f_val = V_field[n] + (M_flow[n])**beta

        for _ in range(NEWTON_MAX_ITER):
            total_newton += 1
            grad_u = Grad @ u_curr
            
            F = A_diff @ u_curr + dt * discrete_Hamiltonian(grad_u) - u_next - dt * f_val
            
            if np.linalg.norm(F, np.inf) < NEWTON_TOL: 
                break
                
            J = A_diff + dt * sp.spdiags(dH_dp(grad_u), 0, Nx, Nx) @ Grad
            u_curr -= spla.spsolve(J, F)
        
        U[n] = u_curr
    return U, total_newton

def solve_fp_forward(U_flow):
    M = np.zeros((Nt + 1, Nx))
    M[0] = m0
    
    for n in range(Nt):
        m_curr = M[n]
        v = - (Grad @ U_flow[n]) 
        
        # Implicit Euler for continuity equation
        Adv_matrix = Grad @ sp.spdiags(v, 0, Nx, Nx)
        A = sp.eye(Nx) - dt*nu*D2 + dt*Adv_matrix
        
        M[n+1] = spla.spsolve(A, m_curr)
        
        # Mass Correction
        M[n+1] = np.maximum(M[n+1], 0)
        M[n+1] /= (np.sum(M[n+1])*Dx + 1e-12)
        
    return M

# ==========================================
# SECTION 5: DAMPED PICARD ITERATION
# ==========================================

start_time_all = time.time()

with open(history_filename, "w") as f_log:
    header = (
        f"{'='*118}\n"
        f"   TRAFFIC LIGHT MFG EXAMPLE (NEWTON TYPE 2: NEWTON-PICARD)\n"
        f"{'='*118}\n"
        f"Parameters:\n"
        f"  xmin = {xmin}, xmax = {xmax}\n"
        f"  T    = {T}, Nx = {Nx}, Nt = {Nt}\n"
        f"  nu   = {nu}, beta = {beta}, penalty = {LIGHT_PENALTY}\n"
        f"Solver Parameters:\n"
        f"  PICARD_MAX_ITER = {PICARD_MAX_ITER}, PICARD_TOL = {PICARD_TOL}\n"
        f"  DAMPING = {DAMPING}\n"
        f"  NEWTON_MAX_ITER = {NEWTON_MAX_ITER}, NEWTON_TOL = {NEWTON_TOL}\n"
        f"Grid Info:\n"
        f"  dt = {dt:.6f}, dx = {Dx:.6f}\n"
        f"{'-'*118}\n"
        f"{'Iter':<5} | {'Abs Err U':<12} | {'Rel Err U':<12} | {'Abs Err M':<12} | {'Rel Err M':<12} | {'Newton It':<10} | {'Time (s)':<10}\n"
        f"{'-'*118}\n"
    )
    print(header, end='')
    f_log.write(header)

    M_flow = np.zeros((Nt + 1, Nx))
    for n in range(Nt+1): M_flow[n] = m0 
    U_flow = np.zeros((Nt + 1, Nx))

    for k in range(1, PICARD_MAX_ITER + 1):
        t0 = time.time()
        
        U_candidate, n_iters = solve_hjb_backward(M_flow)
        U_flow_new = DAMPING * U_candidate + (1 - DAMPING) * U_flow
        M_candidate = solve_fp_forward(U_flow_new)
        M_flow_new = DAMPING * M_candidate + (1 - DAMPING) * M_flow
        
        abs_err_u = np.linalg.norm(U_flow_new - U_flow, 2) * norm_const
        rel_err_u = abs_err_u / (np.linalg.norm(U_flow_new, 2) + 1e-12)
        abs_err_m = np.linalg.norm(M_flow_new - M_flow, 2) * norm_const
        rel_err_m = abs_err_m / (np.linalg.norm(M_flow_new, 2) + 1e-12)
        
        iter_time = time.time() - t0
        
        log_str = f"{k:<5} | {abs_err_u:.4e}   | {rel_err_u:.4e}   | {abs_err_m:.4e}   | {rel_err_m:.4e}   | {n_iters:<10} | {iter_time:.4f}"

        print(log_str)
        f_log.write(log_str + "\n")

        U_flow, M_flow = U_flow_new, M_flow_new
        if rel_err_u < PICARD_TOL and rel_err_m < PICARD_TOL:
            conv_msg = f"{'-'*118}\nCONVERGED at nu={nu} in {k} iterations.\n"
            print(conv_msg)
            f_log.write(conv_msg)
            break

    total_time = time.time() - start_time_all
    time_msg = f"Total Execution Time: {total_time:.4f} seconds.\n"
    print(time_msg)
    f_log.write(time_msg)

# ==========================================
# SECTION 6: VISUALIZATION
# ==========================================

print("Generating Static Plots...")

# --- ONLY CHANGE: INCREASE FONT SIZES FOR PAPER ---
plt.rcParams.update({
    'font.size': 25,            # General text
    'axes.labelsize': 20,       # x and y labels
    'axes.titlesize': 22,       # titles
    'xtick.labelsize': 16,      # x-ticks
    'ytick.labelsize': 16,      # y-ticks
    'legend.fontsize': 16,      # legend
})
# ---

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
target_times = [0.0, T/5, 2*T/5, 3*T/5, 4*T/5, T]
colors = plt.cm.jet(np.linspace(0, 1, len(target_times)))

im = axes[0].imshow(M_flow, aspect='auto', extent=[xmin, xmax, T, 0], cmap='viridis')
axes[0].set_title('Density Evolution (m)')
axes[0].set_xlabel('x')
axes[0].set_ylabel('t')
fig.colorbar(im, ax=axes[0])

for i, t_val in enumerate(target_times):
    n_idx = min(int(round(t_val / dt)), Nt)
    axes[1].plot(x, M_flow[n_idx], color=colors[i], linewidth=2, label=f't={t_val}')
axes[1].set_title('Density Snapshots')
axes[1].legend()
axes[1].grid(True)

for i, t_val in enumerate(target_times):
    n_idx = min(int(round(t_val / dt)), Nt)
    axes[2].plot(x, U_flow[n_idx], color=colors[i], linewidth=2, label=f't={t_val}')
axes[2].set_title('Value Function Snapshots (u)')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.savefig(plot_snapshot_filename)
plt.close()
print(f"Static plots saved to {plot_snapshot_filename}")

# --- Image Sequence Saving ---
print(f"Saving individual frames to {IMG_DIR}/ directory...")

fig_frame, ax_frame = plt.subplots(figsize=(10, 6))

for n in range(Nt + 1):
    ax_frame.clear()
    ax_frame.set_xlim(xmin, xmax)
    ax_frame.set_ylim(0, 10.0)
    ax_frame.set_xlabel("Position (0=Start, 1=End)")
    ax_frame.set_ylabel("Density of Crowd")
    
    t_val = n * dt
    
    if LIGHT_START <= t_val <= LIGHT_END:
        color = 'red'
        alpha = 0.3
        status = "STATUS: RED LIGHT"
    else:
        color = 'green'
        alpha = 0.1
        status = "STATUS: GREEN LIGHT"
        
    ax_frame.add_patch(plt.Rectangle((LIGHT_POS-0.02, 0), 0.04, 10.0, color=color, alpha=alpha))
    
    y_val = M_flow[n]
    ax_frame.plot(x, y_val, 'b-', lw=2)
    ax_frame.fill_between(x, 0, y_val, color='blue', alpha=0.2)
    
    # --- MOVED TIME AND STATUS TO TITLE ---
    ax_frame.set_title(f"Time = {t_val:.2f}s | {status}", color=color, fontweight='bold')
    
    frame_path = os.path.join(IMG_DIR, f"frame_{n:03d}.png")
    plt.savefig(frame_path)

plt.close(fig_frame)
print(f"Successfully saved {Nt + 1} frames to the {IMG_DIR}/ directory.")

# --- Animation Saving ---
print("Generating Animation...")

fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(xmin, xmax)
ax.set_ylim(0, 10.0) 
ax.set_xlabel("Position (0=Start, 1=End)")
ax.set_ylabel("Density of Crowd")

line, = ax.plot([], [], 'b-', lw=2, label='Car Density')
light_patch = plt.Rectangle((LIGHT_POS-0.02, 0), 0.04, 10.0, color='green', alpha=0.1)
ax.add_patch(light_patch)

def init():
    line.set_data([], [])
    ax.set_title("")
    return line, light_patch

def update(frame):
    y = M_flow[frame]
    line.set_data(x, y)
    
    ax.collections.clear()
    ax.add_patch(light_patch) 
    ax.fill_between(x, 0, y, color='blue', alpha=0.2)
    
    t_curr = frame * dt
    
    if LIGHT_START <= t_curr <= LIGHT_END:
        light_color = 'red'
        light_patch.set_color('red')
        light_patch.set_alpha(0.3)
        status_str = "STATUS: RED LIGHT"
    else:
        light_color = 'green'
        light_patch.set_color('green')
        light_patch.set_alpha(0.1)
        status_str = "STATUS: GREEN LIGHT"
        
    # --- MOVED TIME AND STATUS TO TITLE ---
    ax.set_title(f"Time = {t_curr:.2f}s | {status_str}", color=light_color, fontweight='bold')
        
    return line, light_patch

ani = animation.FuncAnimation(fig, update, frames=Nt+1, init_func=init, interval=50, blit=False)
ani.save(video_filename, writer='ffmpeg', fps=30)
print(f"Animation saved to {video_filename}")