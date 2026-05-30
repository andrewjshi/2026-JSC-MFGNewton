Email me at andrewshi@math.berkeley.edu for bugs, questions, requests, support, etc. 

# Code Description

This repository contains the code for the numerical examples for the paper "Newton-Based Solvers for the Finite-Difference Discretization of Mean Field Games" by Mathieu Lauriere and Andrew Shi, 2026. (ADD ARXIV LINK WHEN AVAILABLE)

There are six different examples in this paper: {discounted, ergodic, doublewell, congestionlocal, congestionnonlocal, trafficlight} which could be implemented in dimension {1d, 2d} and with Newton Type 1: Global Newton or Newton Type 2: Newton-Picard {newtontype1, newtontype2}. 

So in principle could be $6x2x2 = 24$ codes with the naming convention example-dimension-newtontype.py (e.g. discounted-1d-newtontype1.py), but certain dimensions or methods are not implemented for all problems. 

| Problem             | 1D Type 1 | 1D Type 2 | 2D Type 1 | 2D Type 2 |
|---------------------|:---------:|:---------:|:---------:|:---------:|
| Discounted          | ✅        | ✅        | ❌        | ❌        |
| Ergodic             | ❌        | ✅        | ❌        | ❌        |
| Double Well         | ✅        | ✅        | ✅        | ✅        |
| Local Congestion    | ✅        | ✅        | ✅        | ✅        |
| Nonlocal Congestion | ❌        | ✅        | ❌        | ✅        |
| Traffic Light       | ✅        | ✅        | ❌        | ❌        |

# Recreating the Paper Figures and Tables

Each table below maps a figure/table in the paper to the script that generates it and the output file that script produces.

> **Note:** every figure in the paper is a TikZ rendering of the underlying numerical output. Running a script reproduces the **Output** listed below — a `.png` snapshot, a `.txt` iteration history (used by TikZ), or a `.csv` for the tables — and the paper's styled figure is generated from that.

### 1. [The Discounted Problem](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-DISCOUNTED)

| Paper | Script | Output |
|---|---|---|
| Figure 2 | [discounted-1d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)) | [type2-discounted_snapshots_T10.0_nu0.5.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)/type2-discounted_snapshots_T10.0_nu0.5.png) |
| Figure 3 | [discounted-1d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)) | [newton-type2-history_discounted1d_T10.0_nu0.5.txt](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)/newton-type2-history_discounted1d_T10.0_nu0.5.txt) |
| Figure 4 | [discounted-1d-newtontype1.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-1%20(Global%20Newton)) | [newton-type1-history_discounted1d_T10.0_nu0.5.txt](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-1%20(Global%20Newton)/newton-type1-history_discounted1d_T10.0_nu0.5.txt) |

*Figure 2 can also be produced by `discounted-1d-newtontype1.py`. For Figures 3 and 4, both the left and right panels are drawn from the linked data file.*

### 2. [The Ergodic Problem](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-ERGODIC)

| Paper | Script | Output |
|---|---|---|
| Figure 6 | [ergodic-1d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-ERGODIC/1D/Newton-Type-2%20(Newton-Picard)) | [ergodic1d_spatial_snapshots_T10.0_nu0.5.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D-ERGODIC/1D/Newton-Type-2%20(Newton-Picard)/ergodic1d_spatial_snapshots_T10.0_nu0.5.png) |
| Figure 7 | [discounted-1d-rho-limit.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/rho-limit) | [squashing_parabolas_T10.0.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/rho-limit/squashing_parabolas_T10.0.png) (left), [loglog_convergence_T10.0.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/rho-limit/loglog_convergence_T10.0.png) (right) |

### 3. [The Double-Well Problem](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL)

| Paper | Script | Output |
|---|---|---|
| Figure 8 | [doublewell-1d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Newton-Type-2%20(Newton-Picard)) | [doublewell1d-plot_T1.0_nu0.1.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Newton-Type-2%20(Newton-Picard)/doublewell1d-plot_T1.0_nu0.1.png) |
| Figure 9 | [doublewell1d_spectral_radius.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Spectral-Radius) | [rho_P.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Spectral-Radius/rho_P.png) (left), [rho_P_theta.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Spectral-Radius/rho_P_theta.png) (right) |
| Figure 10 | [doublewell-2d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/Newton-Type-2%20(Newton-Picard)) | [doublewell2d-plot_T1.0_nu0.1.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/Newton-Type-2%20(Newton-Picard)/doublewell2d-plot_T1.0_nu0.1.png) |
| Table 3 | [timingandfillin.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/timing-fillin-tests/timingandfillin.py) | [timingandfillin_results.csv](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/timing-fillin-tests/timingandfillin_results.csv) |
| Table 4 | [sweep_np_schedule.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/NP-only/sweep_np_schedule.py) | [np2d_schedule_sweep.csv](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/NP-only/np2d_schedule_sweep.csv) |
| Table 5 | [sweep_gn_schedule.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/GN-only/sweep_gn_schedule.py) | [gn2d_schedule_sweep.csv](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/GN-only/gn2d_schedule_sweep.csv) |

### 4. [Local Congestion](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-CONGESTION-LOCAL)

| Paper | Script | Output |
|---|---|---|
| Figure 12 | [congestionlocal-1d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/1D/Newton-Type-2%20(Newton-Picard)) | [localcongestion1d-plot_T1.0_nu0.01.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/1D/Newton-Type-2%20(Newton-Picard)/localcongestion1d-plot_T1.0_nu0.01.png) |
| Figure 13 | [congestionlocal-2d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/2D/Newton-Type-2%20(Newton-Picard)) | [localcongestion2d-contour_T1.0_nu0.01.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/2D/Newton-Type-2%20(Newton-Picard)/localcongestion2d-contour_T1.0_nu0.01.png) |

### 5. [Non-Local Congestion](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20)

| Paper | Script | Output |
|---|---|---|
| Figure 14 | [congestionnonlocal-1d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/1D/Newton-Type-2%20(Newton-Picard)) | [nonlocalcongestion1d-plot_T1.0_nu0.01.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/1D/Newton-Type-2%20(Newton-Picard)/nonlocalcongestion1d-plot_T1.0_nu0.01.png) |
| Figure 15 | [congestionnonlocal-2d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/2D/Newton-Type-2%20(Newton-Picard)) | [nonlocalcongestion2d-contour_T1.0_nu0.01.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/2D/Newton-Type-2%20(Newton-Picard)/nonlocalcongestion2d-contour_T1.0_nu0.01.png) |

### 6. [Traffic Light](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT)

| Paper | Script | Output |
|---|---|---|
| Figure 17 | [trafficlight-1d-newtontype2.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-2%20(Newton-Picard)) (run at θ = 0.3 and θ = 0.1) | [newton-type2-history_trafficlight1d_T1.0_nu0.02.txt](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-2%20(Newton-Picard)/newton-type2-history_trafficlight1d_T1.0_nu0.02.txt) (one per θ) |
| Figure 18 | [trafficlight-1d-newtontype1.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-1%20(Global%20Newton)) | [img/](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-2%20(Newton-Picard)/img) frames: frame_000, 029, 050, 069, 080, 100 .png |

*Figure 18 can also be produced by `trafficlight-1d-newtontype2.py`.*
