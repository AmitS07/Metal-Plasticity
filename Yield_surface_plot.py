# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 21:39:08 2023

@author: Amit Singh
"""
### Code to plot Yield surface ###
#%% 0- Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

#%% # Ploting requirment
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.figsize'] = [9, 9] #default: [6.4, 4.8]
#plt.figure[(9,9)]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'DejaVu Sans'
# plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
# plt.rcParams['text.latex.preamble'] = [r'\boldmath']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'

#%% The Equation to plot Yield Surface 
## Yield function as per the Question 3 in assignment 2, Theroy of plasticity
sigma_o = 200 # Value given
sigma = np.arange(-200,201,10), # to interpolate the data
## Calcuation of tau
tau_p = [] # +ve tau value
tau_n = [] # -ve tau value
for sigma_x in sigma:
    tau_xy_p = sqrt((1/3)*(sigma_o**2 - sigma_x**2))
    tau_xy_n = -1*tau_xy_p
    tau_p.append(tau_xy_p)
    tau_n.append(tau_xy_n)     
tau_p = np.array(tau_p)
tau_n = np.array(tau_n)    
plt.plot(sigma,tau_p,color ='darkblue', marker='o',markersize=5,linewidth =3,label=r'Yield surface: $\sigma$ vs $\tau$')
plt.plot(sigma,tau_n,color ='darkblue', marker='o',markersize=5, linewidth =3)
plt.axvline(x=0, color='grey', linestyle='--')
plt.axhline(y=0, color='grey', linestyle='--')
plt.xlim(-250,250)
plt.ylim(-250,250)

## Calcuation of S11 and S12
S11 = []
S12_p = [] # +ve S12 value
S12_n = [] # -ve S12 value
for stress in sigma:
    P = (1/3)*(stress)
    S_11 = stress -P
    S_12_p = sqrt((1/3)*sigma_o**2 - (1/2)*S_11**2)
    S_12_n = -1* S_12_p
    S11.append(S_11) 
    S12_p.append(S_12_p)   
    S12_n.append(S_12_n)   
plt.plot(S11,S12_p,color ='turquoise', marker='s',markersize=5,linewidth =3,label=r'Yield surface: $S11$ vs $S12$')
plt.plot(S11,S12_n,color ='turquoise', marker='s',markersize=5,linewidth =3)
plt.legend(frameon=False)
plt.xlabel(r"$\sigma$, $S11$ (MPa)")
plt.ylabel(r"$\tau$, $S12$ (MPa)")  
plt.xlim(-250,250)
plt.ylim(-250,250)
plt.savefig("Yield_surface_plot.png", dpi=300,bbox_inches='tight',pad_inches=0.25) 
plt.show()

