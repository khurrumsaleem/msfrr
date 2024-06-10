###############################################################################
# "k_eff_loop.py" runs a sequance of OpenMC simulations. Before each          #
# simulation functions "geometry.py" and "materials" are called to build the  #
# recreated model of the EVOL refference MSFR, at a certain material          #
# temperature (Doppler), and density temperature.                             #
#                                                                             #
# OpenMC is run for a total of 1+3x16=49 times:                               #
#                                                                             #
# 16 eigenvalue simulations with T going from 300K-1800K in steps of 100K.    #
#                                                                             #
# 16 with density T going from 300K-1800K in steps of 100K, while material T  #
#    is constant at the benchmark T=700C. This gives the density contribution #
#    to multiplication facto, relative to the benchmark T.                    #
#                                                                             #
# 16 with material T going from 300K-1800K in steps of 100K, while density T  #
#    is constant at the benchmark T=700C. This gives the "Doppler"            #
#    contribution to multiplication facto, relative to the benchmark T,       #
#    including all effects from the material T changes such as Doppler        #
#    broadening of the neutron interaction cross-sections and hardening/      #
#    softening of the neutron spectrum.                                       # 
#                                                                             #
# 1 with both temperatures equal to the benchmark temperature of 700C.        #
#                                                                             #
#                                                                             #
# Data is saved in files "data/k_eff_XXXXX.txt" with columns:                 #
#                                                                             #
# {particles},{batches},{inactive_cycles},{material_temp},{density_temp},     #
#    {k_combined}                                                             #
#                                                                             #
# OBS: Update "filename" manually to avoid overwriting data.                  #
#                                                                             #
#                                                                             #
# By Morten Nygaard for the DTU student project "Impact of temperature        #
# feedback on reactivity parameters in the molten salt fast reactor"          #
# under the DTU Centre of Nuclear energy technology (ongoing per June 2024).  #
###############################################################################


#%% Initialisation
import numpy as np
import sys
sys.path.append('/opt/openmc/')
import openmc
import openmc.data.multipole
import materials
import geometry

# To visualize computation time, time since start is printed after each simulation.
from time import time
start_time = time()
from datetime import datetime
current_time0 = datetime.now().strftime('%H:%M:%S')
print('Start Time:', current_time0)

#%% Storing data

# OBS: Edit result_filename immediately AFTER starting a sim, to prevent data overwrite
result_filename = 'k_eff_loop_10014.txt'

# OBS: Edit path for saving data
path = "/home/mortennygaard/.cache/yay/openmc-git/msfr/data/"

# OBS: Edit comment 
comment = "multipole=true. create_delayed_neutrons=False. Density from table 1. Inventory from table 7 \n"

# Create or overwrite filename - Comment out to append existing file
with open(path+result_filename, 'w') as f:
   f.write('particles,batches,inactive_cycles,material_temperature_K,density_temperature,k_combined\n')

#%% settings
settings = openmc.Settings()
# minimum of ~100k particles to resolve behaviour of alpha 
# 1M particles would propably be excellent, 300k was used due to limited computation 
particles = 300000
batches = 100
inactive_cycles = 20 # Shannon entropy shown to converge around 10 batches

settings.particles = particles
settings.batches = batches
settings.inactive = inactive_cycles

# Set wether or not to include delayed neutrons, for calculating the effective delayed neutron fraction, beta_eff
#settings.create_delayed_neutrons=False

# The standard temperature interpolation method yielded results appearing unphysical, for unknown reasons. This line sets the interpolation method to "Windowed Multipole", appearing to solve the issue. See OpenMC documentaion.
settings.temperature = {'method':'interpolation','multipole':True}

# Fuel salt outer dimensions, for defining the source and mesh for Shanon Entropy 
r_fuel_salt_out = 225.5/2 + 50.0 + 20.0 + 23.7
h_fuel_salt_out = 225.5/2

# Define Source region
source_area = openmc.stats.Box([-r_fuel_salt_out, -r_fuel_salt_out, -h_fuel_salt_out],[r_fuel_salt_out, r_fuel_salt_out, h_fuel_salt_out],only_fissionable = True)
settings.source = openmc.Source(space=source_area)

# Create xy mesh to calculate Shannon entropy (Converges in ~10 batches)
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = (-r_fuel_salt_out, -r_fuel_salt_out)
entropy_mesh.upper_right = (r_fuel_salt_out, r_fuel_salt_out)
entropy_mesh.dimension = (10, 10)
settings.entropy_mesh = entropy_mesh

# Export all settings to .xml file
settings.export_to_xml()

#%% ### Define sequence for material and density temperatures for loop ###

# Define temperature range in K. 
T0 = np.array([700.+273.15]) # Benchmark T=700C=973.15K
T1 = np.linspace(300.,1800.,16) # In steps of 100K
T2 = T0*np.ones(16) # Constant at benchmark T

# First 
material_temp  =  np.concatenate((T1,T2,T1,T0))
density_temp   =  np.concatenate((T1,T1,T2,T0))

#%% ################# Loop calling OpenMC ########################
# Loop calling the functions "materials.py" and "geometry.py" for varying temperature sequences, for a total of 3x16+1=49 OpenMC simulations giving 49 k_eff data points.

for i in range(np.size(material_temp)):
    
    # Create the materials.xml with the material- and density temperature
    mats, reflector, fuel_salt, b4c, fertile_salt = \
        materials.materials(material_temp[i],density_temp[i]) 
        
    # Create the geometry.xml with the created materials
    geometry.geometry(reflector, fuel_salt, b4c, fertile_salt)
    
    # Run the simulation. OBS: Change number of threads if more are available
    openmc.run(threads=8)
    
    # Load the statepoint file
    sp = openmc.StatePoint('statepoint.'+str(batches)+'.h5')
    
    # Extract the combined k-effective value
    k_combined = sp.k_combined
    
    # Append k_combined and meta data to .txt file
    with open(path+result_filename, 'a') as f:
        f.write(f'\n{particles}, {batches}, {inactive_cycles},  {material_temp[i]}, {density_temp[i]}, {k_combined}')
        
    # Close the statepoint file
    sp.close()
    
    # Print to prompt
    print(r"Material  temperature (C) = "+str(round(material_temp[i]-273.15,2)))
    print(r"Density temperature (C) = "+str(round(density_temp[i]-273.15,2)))
    print(r"Material  temperature (K) = "+str(material_temp[i]))
    print(r"Density temperature (K) = "+str(density_temp[i]))
    print(r"Particles = "+str(particles))
    print(r"Batches = "+str(batches))
    print(r"Inactive Cycles = "+str(inactive_cycles))
    print(r"filename = "+result_filename)
    print(comment)
    print('Start Time:', current_time0)
    current_time = datetime.now().strftime('%H:%M:%S')
    print('Current Time:', current_time)
    print(r"Minutes passed:" + str(round((time()-start_time)/60,2)))
    # End of loop 
    
# Append comment to .txt file
with open(path+result_filename, 'a') as f:
    f.write('\n'+comment)
    
    
###################
### End of code ###
###################
