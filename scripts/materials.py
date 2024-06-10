
###############################################################################
# The function "materials.py" is called by "k_eff_loop.py" to generate the    #
# materials for the recreated OpenMC model of the EVOL reference MSFR as      #
# described in "Neutronic benchmark of the molten salt fast reactor in the    #
# frame of the EVOL and MARS collaborative projects" by Mariya Brovchenko et  #
# al. ("Benchmark article").                                                  #
# See: https://www.epj-n.org/articles/epjn/pdf/2019/01/epjn180012.pdf         #
#                                                                             #
# "materials.py" takes as input the material temperature set as the           #
# temperature of the OpenMC material (Doppler), as well as the density        #
# temperature, which is used to calculate the density of the material.        #
#                                                                             #
# By Morten Nygaard for the DTU student project "Impact of temperature        #
# feedback on reactivity parameters in the molten salt fast reactor"          #
# under the DTU Centre of Nuclear energy technology (ongoing per June 2024).  #
###############################################################################

import openmc

def materials(material_temp,density_temp):
    
    # Starting inventory in kg, rom Benchmark article Table 7
    Th_start = 38281 
    U233_start = 4838 
    
    # Calculate ratios of starting inventory
    Th_ratio_start = Th_start/(Th_start+U233_start)
    U233_ratio_start = U233_start/(Th_start+U233_start)
    
    # Fuel composittion in mol%
    LiF_percent = 0.775
    ThF4_percent = 0.225
    
    # Salt density formula as given in benchmark article
    salt_density = 4.094-8.82e-4*(density_temp-1008)
    
    #%% ###### Defining the materials of the EVOL reference MSFR ######
    
    # Reflector Ni-based Alloy. From benchmark article, Table 2
    reflector = openmc.Material(name='reflector',temperature = material_temp)
    reflector.add_element('Ni',79.432,percent_type='ao')
    reflector.add_element('W',9.976,percent_type='ao')
    reflector.add_element('Cr',8.014,percent_type='ao')
    reflector.add_element('Mo',0.736,percent_type='ao')
    reflector.add_element('Fe',0.632,percent_type='ao')
    reflector.add_element('Ti',0.295,percent_type='ao')
    reflector.add_element('C',0.294,percent_type='ao')
    reflector.add_element('Mn',0.257,percent_type='ao')
    reflector.add_element('Si',0.252,percent_type='ao')
    reflector.add_element('Al',0.052,percent_type='ao')
    reflector.add_element('B',0.033,percent_type='ao')
    reflector.add_element('P',0.023,percent_type='ao')
    reflector.add_element('S',0.004,percent_type='ao')
    reflector.set_density('g/cm3',10)
    
    # Fuel salt. LiF-ThF4-233UF4. From benchmark article Table 3. Atomic percent
    fuel_salt = openmc.Material(name='fuel_salt', temperature = material_temp)
    fuel_salt.add_element('Li', LiF_percent,percent_type='ao', enrichment=99.995, enrichment_target='Li7')
    fuel_salt.add_element('F',LiF_percent+ThF4_percent*4,percent_type='ao')
    fuel_salt.add_element('Th',ThF4_percent*Th_ratio_start,percent_type='ao')
    fuel_salt.add_nuclide('U233',ThF4_percent*U233_ratio_start,percent_type='ao')
    fuel_salt.set_density('g/cm3',salt_density)
    
    # Boron Carbide B4C Protector 
    # Composition of natural boron: 19.8% of 10B and 80.2% of 11B. 
    # The B4C density is set to 2.52 g/cm3.
    b4c = openmc.Material(name = 'b4c', temperature = material_temp)
    b4c.add_nuclide('B10',4.0*.198,percent_type='ao')
    b4c.add_nuclide('B11',4.0*.802,percent_type='ao')
    b4c.add_element('C',1.0,percent_type='ao')
    b4c.set_density('g/cm3',2.52)
    
    # Fertile salt. LiF-ThF4. From benchmark article, Table 3.
    fertile_salt = openmc.Material(name='fertile_salt', temperature = material_temp)
    fertile_salt.add_element('Li',LiF_percent,percent_type='ao', enrichment=99.995, enrichment_target='Li7')
    fertile_salt.add_element('F',LiF_percent+ThF4_percent*4,percent_type='ao')
    fertile_salt.add_element('Th',ThF4_percent,percent_type='ao')
    fertile_salt.set_density('g/cm3',salt_density)
    
    # Export to xml
    mats = openmc.Materials([reflector,fuel_salt,b4c,fertile_salt])
    mats.export_to_xml()
    
    return mats, reflector, fuel_salt, b4c, fertile_salt

###################
### End of code ###
###################
