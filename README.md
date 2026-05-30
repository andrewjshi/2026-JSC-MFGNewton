Email me at andrewshi@math.berkeley.edu for bugs, questions, requests, support, etc. 

# Code Description

This repository contains the code for the numerical examples for the paper "Newton Methods for Finite-Difference Discretizations of Mean Field Games: A Numerical Study" by Mathieu Lauriere and Andrew Shi, 2026. (ADD ARXIV LINK WHEN AVAILABLE)

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

# Recreating the Paper Figures

This explains how to recreate the figures containing the numerical results. 

**1. [The Discounted Problem](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-DISCOUNTED)**
   
Figure 2 - Run [discounted-1d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)) (or discounted-1d-newtontype1.py). The result is given by [type2-discounted_snapshots_T10.0_nu0.5.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)/type2-discounted_snapshots_T10.0_nu0.5.png), but what's in the paper is a TiKz version of this. 

Figure 3 - Both the left and the right figure are created by TiKz using the data in [newton-type2-history_discounted1d_T10.0_nu0.5.txt](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)/newton-type2-history_discounted1d_T10.0_nu0.5.txt) which is created as a result of running [discounted-1d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-2%20(Newton-Picard)). 

Figure 4 - Both the left and the right figure are created by TiKz using the data in [newton-type1-history_discounted1d_T10.0_nu0.5.txt](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-1%20(Global%20Newton)/newton-type1-history_discounted1d_T10.0_nu0.5.txt) which is created as a result of running [discounted-1d-newtontype1.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/Newton-Type-1%20(Global%20Newton)). 

**2. [The Ergodic Problem](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-ERGODIC)**

Figure 6 - Run [ergodic-1d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-ERGODIC/1D/Newton-Type-2%20(Newton-Picard)). The result is given by [ergodic1d_spatial_snapshots_T10.0_nu0.5.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D-ERGODIC/1D/Newton-Type-2%20(Newton-Picard)/ergodic1d_spatial_snapshots_T10.0_nu0.5.png), but what's in the paper is a TiKz version of this. 

Figure 7 - Run [discounted-1d-rho-limit.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-DISCOUNTED/1D/rho-limit). The image on the left is [squashing_parabolas_T10.0.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/rho-limit/squashing_parabolas_T10.0.png) and the image on the right is [loglog_convergence_T10.0.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D-DISCOUNTED/1D/rho-limit/loglog_convergence_T10.0.png), but what's in the paper is a TiKz version of this. 

**3. [The Double-Well Problem](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL)**

Figure 8 - Run [doublewell-1d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Newton-Type-2%20(Newton-Picard)). The result is given by [doublewell1d-plot_T1.0_nu0.1.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Newton-Type-2%20(Newton-Picard)/doublewell1d-plot_T1.0_nu0.1.png). 

Figure 9 - Run [doublewell1d_spectral_radius.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Spectral-Radius). The result is given by [rho_P.png]([https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Newton-Type-2%20(Newton-Picard)/doublewell1d-plot_T1.0_nu0.1.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Spectral-Radius/rho_P.png)) and [rho_P_theta.png](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/1D/Spectral-Radius/rho_P_theta.png), but what's in the paper is a TiKz version of this.
 
Figure 10 - Run [doublewell-2d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/Newton-Type-2%20(Newton-Picard)). The result is given by [doublewell2d-plot_T1.0_nu0.1.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/Newton-Type-2%20(Newton-Picard)/doublewell2d-plot_T1.0_nu0.1.png).

Table 3 - Run [timingandfillin.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/timing-fillin-tests/timingandfillin.py) and the data gets written to [timingandfillin_results.csv](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/timing-fillin-tests/timingandfillin_results.csv). 

Table 4 - Run [sweep_np_schedule.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/NP-only/sweep_np_schedule.py) and the data gets written to [np2d_schedule_sweep.csv](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/NP-only/np2d_schedule_sweep.csv).

Table 5 - Run [sweep_gn_schedule.py](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/GN-only/sweep_gn_schedule.py) and the data gets written to [gn2d_schedule_sweep.csv](https://github.com/andrewjshi/2026-JSC-MFGNewton/blob/main/FINITEDIFF-1D2D-DOUBLEWELL/2D/hybrid-method/GN-only/gn2d_schedule_sweep.csv).

**4. [Local Congestion](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-CONGESTION-LOCAL)**

Figure 12 - Run [congestionlocal-1d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/1D/Newton-Type-2%20(Newton-Picard)). The result is given by [localcongestion1d-plot_T1.0_nu0.01.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/1D/Newton-Type-2%20(Newton-Picard)/localcongestion1d-plot_T1.0_nu0.01.png). 

Figure 13 - [Run congestionlocal-2d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/2D/Newton-Type-2%20(Newton-Picard)). The result is given by [localcongestion2d-contour_T1.0_nu0.01.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D2D-CONGESTION-LOCAL/2D/Newton-Type-2%20(Newton-Picard)/localcongestion2d-contour_T1.0_nu0.01.png). 

**5. [Non-Local Congestion](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20)**

Figure 14 - Run [congestionnonlocal-1d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/1D/Newton-Type-2%20(Newton-Picard)). The result is given by [nonlocalcongestion1d-plot_T1.0_nu0.01.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/1D/Newton-Type-2%20(Newton-Picard)/nonlocalcongestion1d-plot_T1.0_nu0.01.png). 

Figure 15 - Run [congestionnonlocal-2d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/2D/Newton-Type-2%20(Newton-Picard)). The result is given by [nonlocalcongestion2d-contour_T1.0_nu0.01.png](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D2D-CONGESTION-NONLOCAL%20/2D/Newton-Type-2%20(Newton-Picard)/nonlocalcongestion2d-contour_T1.0_nu0.01.png). 

**6. [Traffic Light](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT)**

Figure 17 - Run [trafficlight-1d-newtontype2.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-2%20(Newton-Picard)) with both theta = 0.3 and theta = 0.1. Collect the data from the log file [newton-type2-history_trafficlight_T1.0_nu0.02.txt](https://github.com/andrewjshi/mfg-newton/blob/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-2%20(Newton-Picard)/newton-type2-history_trafficlight1d_T1.0_nu0.02.txt) for each run. What's in the paper is a TiKz version. 

Figure 18 - Run [trafficlight-1d-newtontype1.py](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-1%20(Global%20Newton)) or trafficlight-1d-newtontype2.py. In the [img](https://github.com/andrewjshi/mfg-newton/tree/main/FINITEDIFF-1D-TRAFFICLIGHT/1D/Newton-Type-2%20(Newton-Picard)/img) folder, they are frame_000.png, frame_029.png, frame_050.png, frame_069.png, frame_080.png, frame_100.png. 
