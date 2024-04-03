## Readme File for Metal-Plasticity Repository##
##
Author: Amit Singh, The University of Alabama, Tuscaloosa, USA
Date: 10th Oct 2023
##

This repository contains the basic Python codes and functions to do the following analysis related to Plasticity

1. Code: "eng_true_stress_yield_stress_RO_correction.py" does stress-strain analysis from raw data of load vs displacement,
  (i)- The Eng. stress-strain calculation, preload correction, and plot
  (ii)- The Ramberg-Osgood(RO) correction of the Eng. Stress-strain data and plot
  (iii)- Yield offset correction and yield strength estimation and plot
  (iv)- True stress-strain calculation and plot
Note: The user has to define to five the raw data of "load vs displacement" and yield offset value (0.002, default)

2. Code: "yield_surface_plot.py" does plot the yield surface plot on a given Stress-Strain data

3. Code: "radial_return_alogorithm.py" uses Hencky's equations based on J2 flow theory for given stress states of a given material

4. Code: "incremental_deformation.py" does estimate the incremental stress state at a point for given material and Hardening laws

NOTE: The user should be careful while using these codes. These codes were written for a basic given dataset. These codes do help in understanding the algorithms behind Plasticity theories

