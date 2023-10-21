# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 13:01:33 2023

@author: Amit Singh
"""
#0 Import Libraries
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from math import pi, log


#%% ################### Functions ###################
## 1. Yield function: from Ezra
def yield_point(strain,stress, elastic_modulus, yield_offset=0.002):
    E = elastic_modulus
    stress_y_offset = E*strain -E*yield_offset
    for i in range(len(strain)):
        if stress_y_offset[i]>stress[i]:
            n1 = i
            n0 = n1-1
            break
    # The slope and intercept for offset line
    m1 = E
    c1 = -E*yield_offset
    # The slope and intercept for line form stress-strain
    m2 = (stress[n1]-stress[n0])/(strain[n1]-strain[n0])
    c2 = stress[n0]-(m2*strain[n0])
    # yeild point
    y_strain = (c2-c1)/(m1-m2)
    y_stress = m1*y_strain+c1
    values = [y_strain, y_stress]
    return values
## 2. Nearest point
def nearest_id(data, value):
    array  = np.asarray(data)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

## 3. power law
def power_law(x, K, n):
    return K * np.power(x, n)

#%% ################## Plotting ###################
plt.rcParams.update({'font.size':25})
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.figsize'] = [9,9]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'


#%% ################## Data Analysis ###################
#1 Importing data
#1.1 Defining path
pth = os.getcwd()
pth += '\Assignment_1'
#1.2 importing file or raw data
fle = pd.read_csv(f'{pth}\AL-6061-T6.TXT',skiprows = 4,sep='\t')
#### Load vs Displacment
displacement = fle.iloc[0:480,0] # in inch
load = fle.iloc[0:480,2] # in lb
plt.plot(displacement,load,linewidth = 3,color='C3')
plt.xlim(0,0.10)
plt.ylim(0,3500)
plt.xlabel(r'Displacement ($in)$')
plt.ylabel(r'Load ($lb)$')
plt.legend(['Load'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\load_displacemnt_plot.png')
plt.close()

#%%### Engineering stress strain
initial_guage_length = 1 # in inch
initial_diameter = 0.25 # in inch
initial_area = pi*(initial_diameter/2)**2
eng_strain = displacement/initial_guage_length
eng_stress = (load/initial_area)*(1/1000)
plt.plot(eng_strain, eng_stress,linewidth = 3,color='C2')
plt.xlim(0,0.10)
plt.ylim(0,70)
plt.xlabel(r'Eng. Strain')
plt.ylabel(r'Eng. Stress ($ksi)$')
plt.legend(['Eng. stress'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\Eng_stress_strain_plot.png')
plt.close()
#%%1.4 Elastic Modolus calculation: by fitting linear line in elastic regime
y_upper_E = float(input("Input the upper Stress value for E: "))
y_upper_idx = (eng_stress - y_upper_E).abs().idxmin() # To find out the index, based on given input
x = eng_strain[0:y_upper_idx]
y = eng_stress[0:y_upper_idx]
E,b = np.polyfit(x,y,1)
print(f" Elastic Modolus = {E}")
#%% 1.5 Strain offset correction due to preload
offset_strain = -b/E
eng_strain_corr = eng_strain-offset_strain
plt.plot(eng_strain_corr, eng_stress,linewidth = 3,color='C2')
plt.xlim(0,0.10)
plt.ylim(0,70)
plt.xlabel(r'Eng. Strain')
plt.ylabel(r'Eng. Stress ($ksi)$')
plt.legend(['Eng. stress preload corrected'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\Eng_stress_strain_offset_plot.png')
plt.close()
#%% 1.6 Finding out maximum or UTS and 
UTS = max(eng_stress)
print(f" UTS ={UTS}")
max_idx = np.argmax(eng_stress) # np.argmax, fidn out index of max value
eng_strain_uts = eng_strain_corr[0:max_idx+1]
eng_stress_uts = eng_stress[0:max_idx+1]
plt.plot(eng_strain_uts,eng_stress_uts,linewidth = 3,color='C2')
plt.xlim(0,0.10)
plt.ylim(0,70)
plt.xlabel(r'Eng. Strain')
plt.ylabel(r'Eng. Stress ($ksi)$')
plt.legend(['Eng. stress'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\Eng_stress_strain_UTS_plot.png')
plt.close()
#%% 1.7 Calculating True Stress -strain 
eng_strain_t = np.add(eng_strain_uts,1)
true_strain = np.log(eng_strain_t)
true_stress = eng_stress_uts*(eng_strain_t)
plt.plot(true_strain, true_stress,linewidth = 3,color='blue')
plt.xlim(0,0.10)
plt.ylim(0,70)
plt.xlabel(r'True Strain')
plt.ylabel(r'True Stress ($ksi)$')
plt.legend(['True stress'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\True_stress_strain_plot.png')
plt.close()
# %% ploting both true-eng stress and strain together
#plt.plot(eng_strain_uts,eng_stress_uts, true_strain, true_stress, linewidth = 3)
plt.plot(eng_strain_uts,eng_stress_uts, linewidth = 3,color='C2')
plt.plot(true_strain, true_stress,linewidth = 3,color='blue')
plt.xlim(0,0.10)
plt.ylim(0,70)
plt.xlabel(r'Strain')
plt.ylabel(r'Stress ($ksi)$')
plt.legend(['Eng. stress', 'True stress'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\True_Eng_stress_strain_plot.png')
plt.close()
#%% 1.8 Calculating 0.2% offset
yield_offset = 0.002
yield_stress = yield_point(true_strain,true_stress,E,yield_offset)[1]
print(f" Yield Stress = {yield_stress}")
idx_yield_stress = (true_stress - yield_stress).abs().idxmin()

#% The plotting of yield offset line in true stress-strain plot
offset_y_stress = E*true_strain-E*yield_offset
plt.plot(true_strain, true_stress,linewidth = 3,color='blue')
plt.plot(true_strain,offset_y_stress, linewidth = 3,color='black',linestyle = 'dashed')
plt.xlim(0,0.10)
plt.ylim(0,70)
plt.xlabel(r'True Strain')
plt.ylabel(r'True Stress ($ksi)$')
plt.legend(['True stress', 'Yield offset: 0.2\%'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\True_stress_strain_yield_plot.png')
plt.close()

#%% 1.9 Calculating True plastic strain 
true_plastic_strain = true_strain - (true_stress/E) 
plt.plot(true_plastic_strain, true_stress,linewidth = 3,color='blue')
plt.xlim(0,0.10)
plt.ylim(0,70)
plt.xlabel(r'True Plastic Strain')
plt.ylabel(r'True Stress ($ksi)$')
plt.legend(['True stress'], frameon=False,loc = 'upper right')
plt.savefig(f'{pth}\True_plastic_strain_stress_plot.png')
plt.close()
#%% 1.9 estimation of True plastic strain and stress after yield point and saving
true_plastic_strain_y = true_plastic_strain[idx_yield_stress:-1] # true_plastic_strain_y, this is after yield
true_stress_y = true_stress[idx_yield_stress:-1] # true_stress_y, this is after yield
# Saving data
df = np.column_stack((true_plastic_strain_y, true_stress_y))
np.savetxt('true_stress_srain_plastic_6061.txt', df,delimiter='\t', fmt='%f')
# plt.plot(true_plastic_strain_y, true_stress_y,linewidth = 2,color='blue')
# plt.xlabel(r'True Plastic Strain')
# plt.ylabel(r'True Stress ($ksi)$')
# plt.legend(['True Stress'], frameon=False)
# plt.show()
#%%1.10 curve-fitting and Coefficent of determination (R^2, regression analysis)
# Curve fitting
fit_para, fit_covariance = curve_fit(power_law, true_plastic_strain_y, true_stress_y)
K, n = fit_para
print(f" K = {K}, n = {n}")
# Regression coefficient calculation
predicted_stress = power_law(true_plastic_strain_y, K, n)
R2 = r2_score(true_stress_y, predicted_stress)
print(f"R^2: {R2}")
#%% plotting fitted curve
# true_strain_estimated = true_stress/E+(true_stress/K)^(1/n)
# plt.plot(true_plastic_strain_y, true_stress_y, linewidth = 2,color='C0')
# plt.plot(true_plastic_strain_y,predicted_stress, linewidth = 2,color='C1')
# checking sensitivity of K, n on stress-strain
K = 58.3
n = 0.04
estimated_stress = power_law(true_strain, K, n)
R2 = r2_score(true_stress, estimated_stress)
print(f"R^2: {R2}")

#%% plotting True Stress vs True Plastic Strain and RO fitted

# K_values = [58.3,58.3,58.3,58.3,58.3]
# n_values = [0.02, 0.03,0.04,0.05,0.06]
K_values = [52.3,55.3,58.3,61.3,64.3]
n_values = [0.04, 0.04,0.04,0.04,0.04]
# colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
colors = ['coral', 'lawngreen', 'turquoise', 'mediumslateblue','gray']
# estimated_stress = power_law(true_strain, K, n)
# R2 = r2_score(true_stress, estimated_stress)
# print(f"R^2: {R2}")
plt.plot(true_strain, true_stress,linewidth = 3,color='black')
for K, n, rang in zip(K_values, n_values, colors):
    #%#---------------------------------------
    true_strain_estimated = true_stress/E+(true_stress/K)**(1/n)
    plt.plot(true_strain_estimated,true_stress, linewidth = 3,color=rang)
    plt.xlim(0,0.10)
    plt.ylim(0,70)
    plt.xlabel(r'True Strain')
    plt.ylabel(r'True Stress ($ksi)$')
legend_labels = [f'R-O fitted\n(K = {K}, n = {n})' for K, n in zip(K_values, n_values)]
plt.legend(['Exp.'] + legend_labels, frameon=False, loc='best', fontsize='19')
#plt.legend(['Exp.', f'R-O fitted \n(K = {K}, n = {n})'], frameon=False, loc='upper right')
plt.savefig(f'{pth}\True_stress_strain_RO_fitted_K1.png')
plt.close()

# plt.show()
#%% plotting True Stress vs True Plastic Strain and RO fitted 
plt.plot(true_plastic_strain, true_stress,linewidth = 3,color='black')
for K, n, rang in zip(K_values, n_values, colors):
    true_strain_plastic_estimated = (true_stress/K)**(1/n)
    plt.plot(true_strain_plastic_estimated,true_stress, linewidth = 3,color=rang)
    plt.xlim(0,0.10)
    plt.ylim(0,70)
    plt.xlabel(r'True Plastic Strain')
    plt.ylabel(r'True Stress ($ksi)$')
legend_labels = [f'R-O fitted\n(K = {K}, n = {n})' for K, n in zip(K_values, n_values)]
plt.legend(['Exp.'] + legend_labels, frameon=False, loc='best', fontsize='19')
#plt.legend(['Exp.', f'R-O fitted \n(K = {K}, n = {n})'], frameon=False, loc='upper right')
plt.savefig(f'{pth}\True_stress_strain_plastic_RO_fitted_K.png')
plt.close()



















