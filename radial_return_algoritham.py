# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 16:04:14 2023

@author: Amit Singh

"""
## Note: This Script can also be used for strain calculation via Hencky's Equations
#%% SECTION 0: Importing required libaries
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt   

#%% SECTION 1: Given state of stress 
#1. Sigma (stress) matrix
print("Enters the State of the Stress")
sigma_xx = float(input("Enter Sigma_xx:"))
tau_xy = float(input("Enter Tau_xy:"))
tau_xz = float(input("Enter Tau_xz:"))
sigma_yy = float(input("Enter Sigma_yy:"))
tau_yz = float(input("Enter Tau_yz:"))
sigma_zz = float(input("Enter Sigma_zz:"))
print("Enters the Material's Parameters")
sigma_o = float(input("Enter Yield stress:"))
E = float(input("Enter Elastic Modulus:"))
nu = float(input("Enter Poisson's Ratio':"))
K = float(input("Enter R-O's K value:"))
n = float(input("Enter R-O's n value:"))

#%% SECTION 2: Calculation of Von-Mises stress and factor for radial return"beta"
#2.1: Calculating Shear modulus
G = E/(2*(1+nu))
#2.2 The current(given) stress state
sigma_current = np.array([[sigma_xx,tau_xy,tau_xz],
                          [tau_xy,sigma_yy,tau_yz],
                          [tau_xz,tau_yz,sigma_zz]])

#2.3 Von-Mises stress
sigma_vm = 1/sqrt(2)*sqrt((sigma_xx-sigma_yy)**2+(sigma_yy-sigma_zz)**2+(sigma_zz-sigma_xx)**2
             +6*(tau_xy**2+tau_xz**2+tau_yz**2))

#2.4 current devaitoric stress state and Bouble cheking_beta
P = (sigma_xx+sigma_yy+sigma_zz)/3 # mean stress
I = np.identity(3)
S_current = sigma_current - P*I
# Sij_Sij = np.square(S_current)
# beta_S = sqrt(2/3)*sigma_o/sqrt(np.sum(Sij_Sij))
# sigma_yield = sqrt((3/2)*(np.sum(np.square(S_current))))


#2.5 Elastic strain based on Hooke's Law
epsilon_xx_e = (1/E)*(sigma_xx-nu*(sigma_yy+sigma_zz))
epsilon_yy_e = (1/E)*(sigma_yy-nu*(sigma_xx+sigma_zz))
epsilon_zz_e = (1/E)*(sigma_zz-nu*(sigma_xx+sigma_yy))
gamma_xy_e = tau_xy/G # The gamma_xy = tau_xy/G is a Negineering shear strain
gamma_yz_e = tau_yz/G
gamma_xz_e = tau_xz/G

epsilon_current_e_eng = np.array([[epsilon_xx_e,gamma_xy_e,gamma_xz_e],
                           [gamma_xy_e,epsilon_yy_e,gamma_yz_e],
                           [gamma_xz_e,gamma_yz_e,epsilon_zz_e]]) # This is already in Eng. form need not to multiply shear components by 2


#2.6 Calculating plastic strain based on Hencky's Equations
# 1st R-O, we need to estimate plastic strain based on Ramberg-Osgood
epsilon_current_p_RO = (sigma_vm/K)**(1/n)
epsilon_current_p = (3*S_current*epsilon_current_p_RO)/(2*sigma_vm)

#2.6.1 Plastic Eng. Strain by Hencky
epsilon_current_p_eng = np.array([[epsilon_current_p[0,0],2*epsilon_current_p[0,1],2*epsilon_current_p[0,2]],
                        [2*epsilon_current_p[1,0],epsilon_current_p[1,1],2*epsilon_current_p[1,2]],
                        [2*epsilon_current_p[2,0],2*epsilon_current_p[2,1],epsilon_current_p[2,2]]])

#2.7 total strain
epsilon_current_eng = epsilon_current_e_eng+epsilon_current_p_eng
print(epsilon_current_eng)

#2.8 beta factor
beta = sigma_o/sigma_vm


#%% SECTION 3: The return the current state of stress to yield surface
#3.1: sigma = sigma*beta
sigma_yldsrf = sigma_current*beta # sigma_yldsrf state on yield surface

P_yldsrf = (sigma_yldsrf[0,0]+sigma_yldsrf[1,1]+sigma_yldsrf[2,2])/3 # mean stress
I = np.identity(3)
S_yldsrf = sigma_yldsrf-P_yldsrf*I
#3.2: Yiled stress "on" yield surface
sigma_yield = sqrt((3/2)*(np.sum(np.square(S_yldsrf))))

#3.3: Estimating the elstic and platsic strain for "on" yield surface
# Elastic strain based on Hooke's Law
epsilon_yldsrf_xx_e = (1/E)*(sigma_yldsrf[0,0]-nu*(sigma_yldsrf[1,1]+sigma_yldsrf[2,2]))
epsilon_yldsrf_yy_e = (1/E)*(sigma_yldsrf[1,1]-nu*(sigma_yldsrf[0,0]+sigma_yldsrf[2,2]))
epsilon_yldsrf_zz_e = (1/E)*(sigma_yldsrf[2,2]-nu*(sigma_yldsrf[0,0]+sigma_yldsrf[1,1]))
gamma_yldsrf_xy_e = sigma_yldsrf[0,1]/G
gamma_yldsrf_yz_e = sigma_yldsrf[1,2]/G
gamma_yldsrf_xz_e = sigma_yldsrf[0,2]/G

epsilon_yldsrf_e_eng = np.array([[epsilon_yldsrf_xx_e,gamma_yldsrf_xy_e,gamma_yldsrf_xz_e],
                           [gamma_yldsrf_xy_e,epsilon_yldsrf_yy_e,gamma_yldsrf_yz_e],
                           [gamma_yldsrf_xz_e,gamma_yldsrf_yz_e,epsilon_yldsrf_zz_e]])

#3.4: Calculating plastic strain based on Hencky's Equations
# RO do 1st, we need to estimate plastic strain based on Ramberg-Osgood
epsilon_yldsrf_p_RO = (sigma_yield/K)**(1/n)
epsilon_yldsrf_p = (3*S_yldsrf*epsilon_yldsrf_p_RO)/(2*sigma_yield)

#3.5: "epsilon_yldsrf_e_eng" is obtained by multiplying 2 in to xy,xz,yz components
epsilon_yldsrf_p_eng = np.array([[epsilon_yldsrf_p[0,0],2*epsilon_yldsrf_p[0,1],2*epsilon_yldsrf_p[0,2]],
                        [2*epsilon_yldsrf_p[1,0],epsilon_yldsrf_p[1,1],2*epsilon_yldsrf_p[1,2]],
                        [2*epsilon_yldsrf_p[2,0],2*epsilon_yldsrf_p[2,1],epsilon_yldsrf_p[2,2]]])

#3.6: total eng. strain
epsilon_yldsrf_eng = epsilon_yldsrf_e_eng+epsilon_yldsrf_p_eng


#%% SECTION 4: Incremntal Deformation plasticity: 
#4.1: The remianing or difference between the stress = sigma_current-sigma_yldsrf
delta_sigma = sigma_current-sigma_yldsrf
# Incremental stpes
n_stpes = 1000
d_sigma = delta_sigma/n_stpes
#4.2: evaluation incremental plastic and elastic strain

sigma_yld = []
epsilon_i_yld = []
epsilon_i_yld_eng = []
epsilon_yld_eng = []
for i in range(1,n_stpes+1):
    sigma_yldsrf = sigma_yldsrf+d_sigma
    sigma_yld.append(sigma_yldsrf)
    P = (sigma_yldsrf[0,0]+sigma_yldsrf[1,1]+sigma_yldsrf[2,2])/3 # mean stress
    I = np.identity(3)
    S_yldsrf = sigma_yldsrf-P*I
    sigma_yield = sqrt((3/2)*(np.sum(np.square(S_yldsrf))))
    
    ## incremnetal Plastic strain
    H = K*n*(sigma_yield/K)**((n-1)/n)
    term1 = np.sum(np.multiply(S_yldsrf,d_sigma))
    term2 = ((3/2)*(1/H))/((2/3)*sigma_yield**2)
    factor = term2*term1
    epsilon_i_yldsrf_p = factor*S_yldsrf
    
    ## incremnetal Engineering Elastic strain
    epsilon_i_xx_e = (1/E)*(d_sigma[0,0]-nu*(d_sigma[1,1]+d_sigma[2,2]))
    epsilon_i_yy_e = (1/E)*(d_sigma[1,1]-nu*(d_sigma[0,0]+d_sigma[2,2]))
    epsilon_i_zz_e = (1/E)*(d_sigma[2,2]-nu*(d_sigma[0,0]+d_sigma[1,1]))
    gamma_i_xy_e = d_sigma[0,1]/G
    gamma_i_yz_e = d_sigma[1,2]/G
    gamma_i_xz_e = d_sigma[0,2]/G
    
    epsilon_i_yldsrf_e_eng = np.array([[epsilon_i_xx_e,gamma_i_xy_e,gamma_i_xz_e],
                               [gamma_i_xy_e,epsilon_i_yy_e,gamma_i_yz_e],
                               [gamma_i_xz_e,gamma_i_yz_e,epsilon_i_zz_e]])
    
    ## incremnetal Engineering Plastic strain
    epsilon_i_yldsrf_p_eng = np.array([[epsilon_i_yldsrf_p[0,0],2*epsilon_i_yldsrf_p[0,1],2*epsilon_i_yldsrf_p[0,2]],
                            [2*epsilon_i_yldsrf_p[1,0],epsilon_i_yldsrf_p[1,1],2*epsilon_i_yldsrf_p[1,2]],
                            [2*epsilon_i_yldsrf_p[2,0],2*epsilon_i_yldsrf_p[2,1],epsilon_i_yldsrf_p[2,2]]])
    
    epsilon_i_yldsrf_eng = epsilon_i_yldsrf_e_eng + epsilon_i_yldsrf_p_eng
    epsilon_i_yld_eng.append(epsilon_i_yldsrf_eng)
    
    # this to to add the elastic eng. strain iteratively
    if i ==1:
        epsilon_yldsrf_eng = epsilon_yldsrf_eng + epsilon_i_yldsrf_eng
    else: 
        epsilon_yldsrf_eng = epsilon_yldsrf_eng +epsilon_i_yldsrf_eng
    epsilon_yld_eng.append(epsilon_yldsrf_eng) # This is sum to the elastic strain of stating yield surface to the elastic strain increment
print("The final strain",epsilon_yld_eng[-1])    
## NOTE: To run for different number of increments, before running 
# this section 4,one has to run section 3 just to make sure you are doing 
# the increments for the starting yield condition 
#%% Save data
strain_1000 = epsilon_yld_eng
stress_1000 = sigma_yld
    
#%% 6. Ploting Sigma vs Epsilon 
# ################## Plotting ingredients ###################
plt.rcParams.update({'font.size':25})
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.figsize'] = [9,9]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'

# sixe of the matrix
n_r, n_c = sigma_current.shape # n_r = # rows, n_c = # column
# for j in range(len(colors)):
for nr in range(n_r):
    for nc in range(n_c):
        plt.figure()
        # # for n_steps = 1000
        strain = [epsilon[nr,nc] for epsilon in strain_1000]
        stress = [sigma[nr,nc] for sigma in stress_1000]
        plt.plot(strain, stress,'o-',linewidth = 3,color='cornflowerblue',label ='1000 increments')
        # for n_steps = 100
        strain = [epsilon[nr,nc] for epsilon in strain_100]
        stress = [sigma[nr,nc] for sigma in stress_100]
        plt.plot(strain, stress,'v-',linewidth = 3,color='palegreen',label ='100 increments')
        # # for n_steps = 10
        strain = [epsilon[nr,nc] for epsilon in strain_10]
        stress = [sigma[nr,nc] for sigma in stress_10]
        plt.plot(strain, stress,'s-',linewidth = 3,color='lightcoral',label ='10 increments')
        lb = f'{nr+1}{nc+1}'
        plt.xlabel(rf'$\varepsilon_{{{lb}}}$')
        plt.ylabel(rf'$\sigma_{{{lb}}}$ (MPa)')
        # plt.xlim(0,0.002)
        # plt.ylim(0,300) 
        plt.legend(frameon = False)
        
        plt.savefig(f'D:/Data_E_Drive/University_of_Alabama//Courses at UA//AEM 691//Assignments/midterm/stress_eng_strain_{nr+1}{nc+1}_{n}_new.png', bbox_inches ='tight',pad_inches=0.1)
        # plt.show()
        
plt.close()










    
    

