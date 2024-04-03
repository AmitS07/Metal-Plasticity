# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:49:33 2023

@author: Amit Singh
"""
## This is for the 2nd question of MIDTERM ##

#%% SECTION 0: Importing required libaries
import numpy as np
import pandas as pd
from math import sqrt
import matplotlib.pyplot as plt  
import matplotlib as mpl
import os

#%% Changing directoy
os.chdir('D:/Data_E_Drive/University_of_Alabama/Courses at UA/AEM 691/Assignments/midterm')

#%% SECTION 2: Importing Given Materials parameters
print("Enters the Material's Parameters")
sigma_o = float(input("Enter Yield stress:"))
E = float(input("Enter Elastic Modulus:"))
nu = float(input("Enter Poisson's Ratio':"))
K = float(input("Enter R-O's K value:"))
n = float(input("Enter R-O's n value:"))  #n = 0.4
#2.1: Calculating Shear modulus
G = E/(2*(1+nu))

#%% SECTION1: Importing data
# Note: The raw data was modified for the 2nd and 3rd segment
raw_data = pd.read_excel("midterm_exam.xlsx",sheet_name="Sheet1") 
# #%# 1st segment data: Loading
data_load = raw_data.iloc[1:34,1:3]
# # #%# 2nd segment data: Unloading
# data_load = raw_data.iloc[35:51,1:3]
# # #%# 3rd segment data: Re-laoding
# data_load = raw_data.iloc[52:120,1:3]

# State of stress at yield surface during laoding
#% 1st segment data: For Loading
data_load_yldsrf = data_load.iloc[17] # The yield point at the surface is a point in data corrspond to 150 VM stress
# #%# 2nd segment data: For Unloading/Relaoding
# data_load_yldsrf = data_load.iloc[0] 
# # #%# 3rd segment data: For Relaoding
# data_load_yldsrf = data_load.iloc[0] 


#%%###############################################################################
sigma_xx = data_load_yldsrf[0]
sigma_yy = 0
sigma_zz = 0
tau_yz = 0
tau_xz = 0
tau_xy = data_load_yldsrf[1]

#2.2 The current(given) stress state
sigma_yldsrf = np.array([[sigma_xx,tau_xy,tau_xz],
                          [tau_xy,sigma_yy,tau_yz],
                          [tau_xz,tau_yz,sigma_zz]])

#2.3 Von-Mises stress
sigma_yldsrf_VM = 1/sqrt(2)*sqrt((sigma_xx-sigma_yy)**2+(sigma_yy-sigma_zz)**2+(sigma_zz-sigma_xx)**2
             +6*(tau_xy**2+tau_xz**2+tau_yz**2))

#2.4 current devaitoric stress state and Bouble cheking_beta
P = (sigma_xx+sigma_yy+sigma_zz)/3 # mean stress
I = np.identity(3)
S_yldsrf = sigma_yldsrf - P*I

#2.5 Elastic strain based on Hooke's Law
epsilon_xx_e = (1/E)*(sigma_xx-nu*(sigma_yy+sigma_zz))
epsilon_yy_e = (1/E)*(sigma_yy-nu*(sigma_xx+sigma_zz))
epsilon_zz_e = (1/E)*(sigma_zz-nu*(sigma_xx+sigma_yy))
gamma_xy_e = tau_xy/G
gamma_yz_e = tau_yz/G
gamma_xz_e = tau_xz/G
#Note: Here shear strain commponemet were calculated as of gamma_xy = 2*epsilon_xy
# this equation is by defalut Eng. elastic Strain
epsilon_yldsrf_e_eng = np.array([[epsilon_xx_e,gamma_xy_e,gamma_xz_e],
                           [gamma_xy_e,epsilon_yy_e,gamma_yz_e],
                           [gamma_xz_e,gamma_yz_e,epsilon_zz_e]])  


#2.6 Calculating plastic strain based on Hencky's Equations
# 1st R-O, we need to estimate plastic strain based on Ramberg-Osgood
epsilon_yldsrf_p_RO = (sigma_yldsrf_VM/K)**(1/n)
epsilon_yldsrf_p = (3*S_yldsrf*epsilon_yldsrf_p_RO)/(2*sigma_yldsrf_VM)
epsilon_yldsrf_p_eng = np.array([[epsilon_yldsrf_p[0,0],2*epsilon_yldsrf_p[0,1],2*epsilon_yldsrf_p[0,2]],
                        [2*epsilon_yldsrf_p[1,0],epsilon_yldsrf_p[1,1],2*epsilon_yldsrf_p[1,2]],
                        [2*epsilon_yldsrf_p[2,0],2*epsilon_yldsrf_p[2,1],epsilon_yldsrf_p[2,2]]])

#2.7 total engineering strain
epsilon_yldsrf_eng = epsilon_yldsrf_e_eng+epsilon_yldsrf_p_eng

#2.9 beta factor
beta = sigma_o/sigma_yldsrf_VM

#%% Now importing remaning data for incremental loading
# State of stress from yield surface to final laoding
#%# 1st segment data: Loading
# data_load_incr = data_load.iloc[18:,:] # The yield point at the surface is a point in data corrspond to 150 VM stress
# # # %# 2nd segment data: For Unloading
data_load_incr = data_load.iloc[1:,:]
# #%# 3rd segment data: For Re-loading
# data_load_incr = data_load.iloc[1:,:]
#####################################
sigma_load_incr = []
epsilon_i_incr_tt_eng = [] # tt = total
epsilon_load_incr_eng = []
for i in range(len(data_load_incr)):
    data_load_i = data_load_incr.iloc[i]
    sigma_xx = data_load_i[0]
    sigma_yy = 0
    sigma_zz = 0
    tau_yz = 0
    tau_xz = 0
    tau_xy = data_load_i[1]
    #2.2 The current(given) stress state
    sigma_load_incr_i = np.array([[sigma_xx,tau_xy,tau_xz],
                              [tau_xy,sigma_yy,tau_yz],
                              [tau_xz,tau_yz,sigma_zz]])
    sigma_load_incr.append(sigma_load_incr_i)
    

    #2.4 current devaitoric stress state and Bouble cheking_beta
    P_incr = (sigma_xx+sigma_yy+sigma_zz)/3 # mean stress
    S_load_incr = sigma_load_incr_i - P_incr*I
    
    ## Yield 
    sigma_yield_incr = sqrt((3/2)*(np.sum(np.square(S_load_incr))))
    
    ## Now calculating incremntal change in stress state
    # print(sigma_load_incr)
    if i ==0:
        d_sigma_incr = sigma_load_incr_i-sigma_yldsrf
    else:
        d_sigma_incr = sigma_load_incr[i] - sigma_load_incr[i-1]
    # d_sigma_incr.append(d_sigma_i_incr)
    
    ## incremnetal Plastic strain
    #%## This is for the transion laoding path (b/w tension and compression)
    if sigma_yield_incr == 0:
        factor = 0
        eplsilon_i_incr_p = factor*S_load_incr
    else:
        H = K*n*(sigma_yield_incr/K)**((n-1)/n)
        term1 = np.sum(np.multiply(S_load_incr,d_sigma_incr))
        term2 = ((3/2)*(1/H))/((2/3)*sigma_yield_incr**2)
        factor = term2*term1
        epsilon_i_incr_p = factor*S_load_incr
        epsilon_i_incr_p_eng = np.array([[epsilon_i_incr_p[0,0],2*epsilon_i_incr_p[0,1],2*epsilon_i_incr_p[0,2]],
                                [2*epsilon_i_incr_p[1,0],epsilon_i_incr_p[1,1],2*epsilon_i_incr_p[1,2]],
                                [2*epsilon_i_incr_p[2,0],2*epsilon_i_incr_p[2,1],epsilon_i_incr_p[2,2]]])

    
    ## incremnetal Elastic strain
    epsilon_i_xx_e = (1/E)*(d_sigma_incr[0,0]-nu*(d_sigma_incr[1,1]+d_sigma_incr[2,2]))
    epsilon_i_yy_e = (1/E)*(d_sigma_incr[1,1]-nu*(d_sigma_incr[0,0]+d_sigma_incr[2,2]))
    epsilon_i_zz_e = (1/E)*(d_sigma_incr[2,2]-nu*(d_sigma_incr[0,0]+d_sigma_incr[1,1]))
    gamma_i_xy_e = d_sigma_incr[0,1]/G
    gamma_i_yz_e = d_sigma_incr[1,2]/G
    gamma_i_xz_e = d_sigma_incr[0,2]/G
    
    epsilon_i_incr_e_eng = np.array([[epsilon_i_xx_e,gamma_i_xy_e,gamma_i_xz_e],
                               [gamma_i_xy_e,epsilon_i_yy_e,gamma_i_yz_e],
                               [gamma_i_xz_e,gamma_i_yz_e,epsilon_i_zz_e]])
    # Total incremental strain = incremental(Ealstic and Palstic)Strain
    epsilon_i_incr_eng = epsilon_i_incr_e_eng + epsilon_i_incr_p_eng
    epsilon_i_incr_tt_eng.append(epsilon_i_incr_eng)
    
    
    # This to to add the elastic eng. strain iteratively
    ## I HAVE BIG DOUBT IN THIS SECTION, does "elastic_strain_yield_surface" has to be added on each increemnt?
    if i ==0:
        epsilon_incr_eng = epsilon_yldsrf_eng + epsilon_i_incr_eng
    else: 
        epsilon_incr_eng = epsilon_incr_eng +epsilon_i_incr_eng
    
    epsilon_load_incr_eng.append(epsilon_incr_eng)
    
#%%# Saving variable 
# Stress
np.savez('reloading_stress.npz', *sigma_load_incr)
stress_11_12_reload = sigma_load_incr
# Strain
np.savez('reloading_strain.npz', *epsilon_load_incr_eng)

#%% Merging data for stress
loading_stress = np.load('loading_stress.npz')
unloading_stress = np.load('unloading_stress.npz')
reloading_stress = np.load('reloading_stress.npz')
# Merging data for strain
loading_stress = np.load('loading_stress.npz')
unloading_stress = np.load('unloading_stress.npz')
reloading_stress = np.load('reloading_stress.npz')

# Extract matrices from loaded data sets
stress_1 = [loading_stress[key] for key in loading_stress.files]
stress_2 = [unloading_stress[key] for key in unloading_stress.files]
stress_3 = [reloading_stress[key] for key in reloading_stress.files]

# Merging data for strain
loading_strain = np.load('loading_strain.npz')
unloading_strain = np.load('unloading_strain.npz')
reloading_strain = np.load('reloading_strain.npz')
# Merging data for strain
loading_strain = np.load('loading_strain.npz')
unloading_strain = np.load('unloading_strain.npz')
reloading_strain = np.load('reloading_strain.npz')

# Extract matrices from loaded data sets
strain_1 = [loading_strain[key] for key in loading_strain.files]
strain_2 = [unloading_strain[key] for key in unloading_strain.files]
strain_3 = [reloading_strain[key] for key in reloading_strain.files]

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
n_r, n_c = sigma_yldsrf.shape # n_r = # rows, n_c = # column
# for j in range(len(colors)):
for nr in range(n_r):
    for nc in range(n_c):
        plt.figure()
        # plt.xlim(-0.05,0.05)
        # plt.ylim(-300,300)
        # strain = [epsilon[nr,nc] for epsilon in epsilon_load_incr_eng]
        strain1 = [epsilon[nr,nc] for epsilon in strain_1]
        # print(strain)
        # stress = [sigma[nr,nc] for sigma in sigma_load_incr]
        stress1 = [sigma[nr,nc] for sigma in stress_1]
        plt.plot(strain1, stress1, 'o-',linewidth = 3,color='cornflowerblue',label ='Loading:Tension-Shear')
        ##### Color map for laod and plot
        # color_scale1 = stress1 
        # norm = mpl.colors.Normalize(vmin = min(color_scale1), vmax = max(color_scale1))
        # cmap = mpl.cm.Blues
        # plt.scatter(strain1, stress1,color=cmap(norm(color_scale1)), label ='Loading')
        #####
        lb = f'{nr+1}{nc+1}'
        plt.xlabel(rf'$\varepsilon_{{{lb}}}$')
        plt.ylabel(rf'$\sigma_{{{lb}}}$ (MPa)')
        ## Unload
        strain2 = [epsilon[nr,nc] for epsilon in strain_2]
        stress2 = [sigma[nr,nc] for sigma in stress_2]
        plt.plot(strain2, stress2, 's-',linewidth = 3,color='lawngreen',label ='Unloading:Shear')
        ##### Color map for unlaod
        # color_scale1 = stress2 
        # norm = mpl.colors.Normalize(vmin = min(color_scale1), vmax = max(color_scale1))
        # cmap = mpl.cm.Greens
        # plt.scatter(strain2, stress2,color=cmap(norm(color_scale1)),label ='Unloading')
        ###plot
        # Reload
        strain3 = [epsilon[nr,nc] for epsilon in strain_3]
        stress3 = [sigma[nr,nc] for sigma in stress_3]
        plt.plot(strain3, stress3, 'v-',linewidth = 3,color='crimson',label ='Unloading:Tension')
        ##### Color map for relaod and plot
        # color_scale1 = stress3 
        # norm = mpl.colors.Normalize(vmin = min(color_scale1), vmax = max(color_scale1))
        # cmap = mpl.cm.Reds
        # plt.scatter(strain3, stress3,color=cmap(norm(color_scale1)),label ='Reloading')
        ######
        plt.legend(frameon = False)
        plt.grid(True)
        # plt.show()
        plt.savefig(f'D:/Data_E_Drive/University_of_Alabama//Courses at UA//AEM 691//Assignments/midterm/Prbl2_stress_eng_strain_{nr+1}{nc+1}_{n}.png',bbox_inches='tight',pad_inches=0.1)        
        plt.close()    
    

#%% To plot loading segments
# #%# 1st segment data: Loading
s_x_load1 = raw_data.iloc[1:34,1]
s_xy_load1 = raw_data.iloc[1:34,2]
plt.plot(s_x_load1, s_xy_load1, 'o-',linewidth = 3,color='cornflowerblue',label ='Loading:Tension-Shear')
#%# 2nd segment data: Unloading
s_x_load2 = raw_data.iloc[35:51,1]
s_xy_load2 = raw_data.iloc[35:51,2]
plt.plot(s_x_load2, s_xy_load2, 's-',linewidth = 3,color='lawngreen',label ='Unloading:Shear')
#%# 3rd segment data: Re-laoding
s_x_load3 = raw_data.iloc[52:120,1]
s_xy_load3 = raw_data.iloc[52:120,2]
plt.plot(s_x_load3, s_xy_load3, 'v-',linewidth = 3,color='crimson',label ='Unloading:Tension')

plt.plot(117.117, 54.11, '*', markersize=15, color='k', label=' Point on yield surface')


plt.legend(frameon = False)
plt.grid(True)
plt.xlabel(r"$\sigma_{11}$ (MPa)") # here 'r' is for to make string raw
plt.ylabel(r"$\tau_{12}$ (MPa)")
# plt.show()
plt.savefig("D:/Data_E_Drive/University_of_Alabama//Courses at UA//AEM 691//Assignments/midterm/load_history1.png",bbox_inches='tight',pad_inches=0.1)        
plt.close()   

#%% To plot Stress_11 vs Stress_12
stress_11_12 = np.concatenate((stress_11_12_load, stress_11_12_unload, stress_11_12_reload))
legend_flag = False

for j in range(len(stress_11_12)):
    if j < 15:
        sigma_11 = stress_11_12[j][0,0]
        sigma_12 = stress_11_12[j][0,1]
        plt.plot(sigma_11,sigma_12,'o-',linewidth = 3,color='cornflowerblue',label ='Loading:Tension-Shear')
    elif 14<j<30:
        sigma_11 = stress_11_12[j][0,0]
        sigma_12 = stress_11_12[j][0,1]
        plt.plot(sigma_11,sigma_12,'s-',linewidth = 3,color='lawngreen',label ='Unloading:Shear')

    else:

        sigma_11 = stress_11_12[j][0,0]
        sigma_12 = stress_11_12[j][0,1]
        plt.plot(sigma_11,sigma_12,'v-',linewidth = 3,color='crimson',label ='Unloading:Tension')
    # To add legend only once
    
plt.plot(117.117, 54.11, '*', markersize=15, color='k', label=' Point on yield surface')     
plt.grid(True)
plt.xlabel(r"$\sigma_{11}$ (MPa)") # here 'r' is for to make string raw
plt.ylabel(r"$\tau_{12}$ (MPa)")  
plt.savefig("D:/Data_E_Drive/University_of_Alabama//Courses at UA//AEM 691//Assignments/midterm/load_history1_code_{n}.png",bbox_inches='tight',pad_inches=0.1)    
# plt.show()