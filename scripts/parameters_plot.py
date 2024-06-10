###############################################################################
# Gennerates plots from data generated specifically by "k_eff_loop.py".       # 
# Multiplication factor k for ref, Doppler, density and sum is converted to   #
# reactivity, rho, which is numerically differentiated to give the            #
# temperature coefficient of reactivity, alpha (for ref, Doppler, density and #
# sum). For a better understanding of this code, please refer to the          #
# "README.md" file for a project description.                                 #
#                                                                             #
# By Morten Nygaard for the DTU student project "Impact of temperature        #
# feedback on reactivity parameters in the molten salt fast reactor"          #
# under the DTU Centre of Nuclear energy technology (ongoing per June 2024).  #
#                                                                             #
# Calculations of standard error propagation are shown in project report.     #
###############################################################################

import os
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 23})
plt.rcParams['text.usetex'] = True # Requires Latex installed
plt.rcParams['font.family'] = 'serif'

#%% Load text file to DataFrame and calc rho and coefs with error

# OBS - Set working directory for the data, "...\.cache\yay\openmc-git\msfr\data"
os.chdir(r"\\wsl.localhost\Manjaro\home\mortennygaard\.cache\yay\openmc-git\msfr\data")

data_filename = r"k_eff_loop_10013.txt"
df = pd.read_csv(data_filename)

# Split multiplication factor string into floats in seperate columns 
df[['k_combined_value','k_combined_error']] = df['k_combined'].str.split('\\+/-',expand=True).astype(float)
L = np.size(df['k_combined_value']) # number of elements
T0 = 700+273.15 # benchmark temperature

# Calculate rho and standard error propagation
df['rho_combined']=( 1-1/df['k_combined_value'] )*1e5 # in pcm
df['rho_combined_error']=df['k_combined_error']/(df['k_combined_value']**2)*1e5

######## Divide data set into appropriately named collumns ########

# DataFrame contains: [start:end+1]
# Rows 0-15  --> [0:16] has T varying together. (Ref)
# rows 16-31 --> [16:32] has constant material temp 700C, varying density temp (density)
# rows 32-47 --> [32:48] has constant density temp 700C, varying material temp (Doppler)
# Row 48     --> [48] is 700C benchmark. 
# Length = 49

# Rename rho_j, j in D,d,R for Doppler, density, and Ref 
df['rho_Doppler'] = df['rho_combined'][32:48]; 
df['rho_Doppler']=df['rho_Doppler'].shift(-32) # move to first element
df['rho_density'] = df['rho_combined'][16:32]; 
df['rho_density']=df['rho_density'].shift(-16) # move to first element
df['rho_Ref'] = df['rho_combined'][0:16]; 

# Calc rho_s, for sum
df['rho_sum'] = df['rho_Doppler']+df['rho_density']-df['rho_combined'][48]

# Rename error on rho_j, j in D,d,R
df['rho_Doppler_error'] = df['rho_combined_error'][32:48]; df['rho_Doppler_error']=df['rho_Doppler_error'].shift(-32) # move to first element
df['rho_density_error'] = df['rho_combined_error'][16:32]; df['rho_density_error']=df['rho_density_error'].shift(-16) # move to first element
df['rho_Ref_error'] = df['rho_combined_error'][0:16]; 

# Calc standard error propagation on rho_s
df['rho_sum_error'] = np.sqrt(df['rho_Doppler_error']**2+\
                              df['rho_density_error']**2+\
                              df['rho_combined_error'][48]**2)


#%% Plot rho_j, j in D,d,s,R for Doppler, density, sum, and Ref 

# Define correlation functions for curve fitting

def log_fit(x,a,b): # Doppler
    return a + b * np.log10(x)

def polynom_fit_2nd_order(x,a,b,c): # density
    return a * x**2 + b * x + c

def linear_fit(x,a,b): # Ref and sum
    return a * x + b

# Initiate figure
fig, ax = plt.subplots(figsize=(10,7)); 
T = np.linspace(300,1800)

# Initiate fit range from 800K-1800K
fit_from = 5
fit_to = 16
marker_size = 5

################## Doppler ###########

# Curve fit
popt,pcov = curve_fit(log_fit,
                      df['material_temperature_K'][fit_from:fit_to],
                      df['rho_Doppler'][fit_from:fit_to], 
                      sigma=df['rho_Doppler_error'][fit_from:fit_to],
                      absolute_sigma=True)
perr = np.sqrt(np.diag(pcov)) # Weighted error on fit parameters

# Plot data and fit
plt.errorbar(df['material_temperature_K'][0:16],
             df['rho_Doppler'][0:16],
             df['rho_Doppler_error'][0:16],
             fmt='bo',ms=0,ecolor='red',capsize=marker_size)
plt.plot(T,log_fit(T,popt[0],popt[1]), color='red',label='Doppler: log fit')

############### density ##########

# Curve fit
popt,pcov = curve_fit(polynom_fit_2nd_order, 
                      df['material_temperature_K'][fit_from:fit_to],
                      df['rho_density'][fit_from:fit_to],
                      sigma=df['rho_density_error'][fit_from:fit_to],
                      absolute_sigma=True)
perr = np.sqrt(np.diag(pcov)) 

# Plot data and fit
plt.errorbar(df['material_temperature_K'][0:16],
             df['rho_density'][0:16],
             df['rho_density_error'][0:16], 
             fmt='bo',ms=0,ecolor='blue',capsize=marker_size)
plt.plot(T,polynom_fit_2nd_order(T,popt[0],popt[1],popt[2]), 
         color='blue',label=r'density: 2$^{\mathrm{nd}}$ order polynom fit')

########## Sum ########

# Curve fit
popt,pcov = curve_fit(linear_fit, 
                      df['material_temperature_K'][fit_from:fit_to],
                      df['rho_sum'][fit_from:fit_to],
                      sigma=df['rho_sum_error'][fit_from:fit_to],
                      absolute_sigma=True)
perr = np.sqrt(np.diag(pcov)) 

# Plot data and fit
plt.errorbar(df['material_temperature_K'][0:16],
             df['rho_sum'][0:16],
             df['rho_sum_error'][0:16],
             fmt='bo',ms=0,ecolor='purple',capsize=marker_size)
plt.plot(T,linear_fit(T,popt[0],popt[1]), 
         color='purple',label='sum: linear fit')

############## Ref ##########

# No fit
# Plot data
plt.errorbar(df['material_temperature_K'][0:16],
             df['rho_Ref'][0:16], 
             df['rho_Ref_error'][0:16],
             fmt='bo',ms=0,ecolor='green',capsize=marker_size,label='Ref')

# Format plot
ax.legend() 
plt.ylabel(r'$\rho$ [pcm]',rotation=90); 
plt.xlabel(r'T [K]');
plt.title(r'data file: '+data_filename+\
          ',  [batches,inactive_cycles,particles]=['+str(df['batches'][1])+\
              ','+str(df['inactive_cycles'][1])+','+str(df['particles'][1])+']'
              ,fontsize=17)
plt.xlim(300,1800)
plt.ylim(-9000,5000)
plt.show()

#%% Saves image if enabled
# OBS - change directory

#os.chdir(r'C:\Users\s182487\OneDrive - Danmarks Tekniske Universitet\Dokumenter\DTU pc jan-24\Impact of temperature feedback on reactivity in the Molten Salt Fast Reactor\plots_constant')
fig.savefig('rho_plot_'+data_filename[0:-4]+'.png', format='png', dpi=300)

#%% ############# Calc coefs and standard error propagation ############

# Note: T steps for coef, T=[350,450...,1750] are in between T steps for rho, 
# T=[300,400,...,1800], and thus 1 element shorter.

df['T_coef'] = df['material_temperature_K'][0:15]+50

# The coefficients are calculated by numerical differentiation of rho:
# alpha{i} = (rho_{i+1}-rho_{i})/(T_{i+1}-T_{i}) = rho.diff()/100K

# Calc coef Doppler and standard error propagation
df['coef_Doppler'] = df['rho_Doppler'].diff()/100
df['coef_Doppler'] = df['coef_Doppler'].shift(-1) # move to first element
df['coef_Doppler_error'] = np.concatenate((np.array(
        np.sqrt([sum(pair) for pair in zip(
        df['rho_Doppler_error'][0:15]**2,
        df['rho_Doppler_error'][0:15]**2)])/100),
        np.ones(L-15)*np.nan)) # NaN's are for matching length to df

# Calc coef density and standard error propagation
df['coef_density'] = df['rho_density'].diff()/100
df['coef_density'] = df['coef_density'].shift(-1)
df['coef_density_error'] = np.concatenate((np.array(
        np.sqrt([sum(pair) for pair in zip(
        df['rho_density_error'][0:15]**2,
        df['rho_density_error'][0:15]**2)])/100),
        np.ones(L-15)*np.nan)) 

# Calc coef sum and standard error propagation
df['coef_sum'] = df['rho_sum'].diff()/100
df['coef_sum'] = df['coef_sum'].shift(-1)
df['coef_sum_error'] = np.concatenate((np.array(
        np.sqrt([sum(pair) for pair in zip(
        df['rho_sum_error'][0:15]**2,
        df['rho_sum_error'][0:15]**2)])/100),
        np.ones(L-15)*np.nan)) 

# Calc coef Ref and standard error propagation
df['coef_Ref'] = df['rho_Ref'].diff()/100
df['coef_Ref'] = df['coef_Ref'].shift(-1)
df['coef_Ref_error'] = np.concatenate((np.array(
        np.sqrt([sum(pair) for pair in zip(
        df['rho_Ref_error'][0:15]**2,
        df['rho_Ref_error'][0:15]**2)])/100),
        np.ones(L-15)*np.nan))

#%% Fit and plot coefficients

# Curve Fitting

def reciprocal_fit(x,a,b): # Doppler (derrivative of a log is a reciprocal)
    return a / (x * np.log(10)) + b


def linear_fit(x,a,b): # density (derriv. of 2nd ord. polynom., is 1st order)
    return a * x + b

# Constant fit yilds approximate value of alpha_sum, 
# assumed approximately constant from 800K-1800K
def linear_fit_a_0(x,b):
    return 0*x+b

# Initiate plot
fig, ax = plt.subplots(figsize=(10,7)); 

############### Doppler #########

# Curve fit
popt,pcov = curve_fit(
    reciprocal_fit,
    df['T_coef'][fit_from:fit_to-1],
    df['coef_Doppler'][fit_from:fit_to-1],
    sigma=df['coef_Doppler_error'][fit_from:fit_to-1],
    absolute_sigma=True)
perr = np.sqrt(np.diag(pcov)) 

# Plot fit and data
plt.plot(T,reciprocal_fit(T,*popt),
         color='red',label='Doppler: reciprocal fit')
plt.errorbar(df['T_coef'][0:15],df['coef_Doppler'][0:15],df['coef_Doppler_error'][0:15],fmt='bo',ms=0,ecolor='red',capsize=marker_size)
    
############### density #########

# Curve fit
popt,pcov = curve_fit(
    linear_fit,
    df['T_coef'][fit_from:fit_to-1],
    df['coef_density'][fit_from:fit_to-1],
    sigma=df['coef_density_error'][fit_from:fit_to-1],
    absolute_sigma=True)
perr = np.sqrt(np.diag(pcov)) 

# Plot fit and data
plt.plot(T,linear_fit(T,*popt),
         color='blue',label='density: linear fit')
plt.errorbar(df['T_coef'][0:15],df['coef_density'][0:15],df['coef_density_error'][0:15],fmt='bo',ms=0,ecolor='blue',capsize=marker_size)

############### sum #########

# Curve fit
popt,pcov = curve_fit(
    linear_fit_a_0,
    df['T_coef'][fit_from:fit_to-1],
    df['coef_sum'][fit_from:fit_to-1],
    sigma=df['coef_sum_error'][fit_from:fit_to-1],
    absolute_sigma=True)
perr = np.sqrt(np.diag(pcov)) 

# Print fit paramters (b=alpha)
print('\ncoef_sum fit parameters ')
print('[b]:             '+str(popt))
print('[b_error]:       '+str(perr))

# Plot fit and data
plt.plot(T,linear_fit_a_0(T,*popt),
         color='purple',label='sum: linear fit')
plt.errorbar(df['T_coef'][0:15],df['coef_sum'][0:15],df['coef_sum_error'][0:15],fmt='bo',ms=0,ecolor='purple',capsize=marker_size)
    
############### Ref #########

# Curve fit
popt,pcov = curve_fit(
    linear_fit_a_0,
    df['T_coef'][fit_from:fit_to-1],
    df['coef_Ref'][fit_from:fit_to-1],
    sigma=df['coef_Ref_error'][fit_from:fit_to-1],
    absolute_sigma=True)
perr = np.sqrt(np.diag(pcov)) 

# Print fit paramters (b=alpha)
print('\ncoef_Ref fit parameters ')
print('[b]:             '+str(popt))
print('[b_error]:       '+str(perr))

# Plot data, no fit
#plt.plot(T,linear_fit_a_0(T,*popt),color='green')
plt.errorbar(df['T_coef'][0:15],df['coef_Ref'][0:15],df['coef_Ref_error'][0:15],fmt='bo',ms=0,ecolor='green',capsize=marker_size,label='Ref')
    
########## Plot formatting ########

plt.plot(T,T*0,'k') # Plot 0-line to visualize alpha<0
plt.legend()
plt.ylabel(r'$\frac{\mathrm{d\rho}}{\mathrm{dK}}$ [pcm/K]',rotation=90); 
plt.xlabel(r'T [K]');
plt.title(r'data file: '+data_filename+'  [batches,inactive_cycles,particles]=['+str(df['batches'][1])+','\
          +str(df['inactive_cycles'][1])+','+str(df['particles'][1])+']',fontsize=17)
plt.xlim(300,1800)
plt.ylim(-13.99,.5)
plt.show()

#%%  Saves coef figure if enabled - change directory
#os.chdir(r'C:\Users\s182487\OneDrive - Danmarks Tekniske Universitet\Dokumenter\DTU pc jan-24\Impact of temperature feedback on reactivity in the Molten Salt Fast Reactor\plots_constant')
fig.savefig('coef_plot_'+data_filename[0:-4]+'.png', format='png', dpi=300)

#%% Export dataframe if enabled
#df.to_csv(data_filename,sep=',',index=False,encoding='utf-8')


###################
### End of code ###
###################
