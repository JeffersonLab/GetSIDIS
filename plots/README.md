# SIDIS Events were generated based on the following kinematics:
  ## P_ebeam = 10 GeV/c, P_ion = 600 GeV/c for the C12 target
  ## Phase-Space are:
  * Electrons:   
```js              
      Pe_min = 0.5 GeV/c,      Pe_max = 30.0 GeV/c
      Theta_min = 0.0 degrees, Theta_max = 140 degrees
      Phi_min = 0.0 degrees,   Phi_max = 360 degrees
```
  * Hadrons (pin+/- in this case):
```js
      Ph_min = 0.0 GeV/c,      Pe_max = 10.0 GeV/c
      Theta_min = 0.0 degrees, Theta_max = 180 degrees
      Phi_min = 0.0 degrees,   Phi_max = 360 degrees
```

# Event Counts are based on Luminisity = 10^33 / A, where A = 12. Assuming full EIC acceptance and with one day of beam time.

# Plotting:
## Pi+ Acceptance in ./pip_acc (similar for Pi- ):
```js 
    c12_pip_Q2_x_log_A600.pdf:  Q2 vs. x, both axies are in log-scale
    c12_pip_acc_A600.pdf: Momentum vs. Theta for both electrons and pions, respectively.
    c12_pip_Q2_z_A600_xbin.pdf: 2D plot of Q2(log) vs. Z in a small x-bin (0.008<x<0.012). I binned 7 Q2-bins and 7 z-bins and calculate how many counts in each bin with 1-day of EIC running.
```
## Per Christian's request, I made few other plots here:
```js
    c12_diff_Q2_z_A600_xbin.pdf: Simila to the pi+ and pi- individual plots, but here I calculated the differences of pi+ counts to pi- counts.
    c12_pion_multiplicity_z.png: "Multiplicity" vs. Z, but here "Multiplicity" is just simply the ratio of SIDIS XS to inclusive XS. 
                                  Three plots: pi+, pi-, and (pi+)-(pi-)
    c12_pion_pt.png: 1D distributions of pt.
```
    


   
                   
