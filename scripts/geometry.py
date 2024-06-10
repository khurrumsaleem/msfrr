###############################################################################
# "geometry.py" is called by k_eff_loop.py" to build recreated OpenMC model   #
# of the EVOL reference MSFR as described in "Neutronic benchmark of the      #
# molten salt fast reactor in the frame of the EVOL and MARS collaborative    #
# projects" by Mariya Brovchenko et al. ("Banchmark article")                 #
# See: https://www.epj-n.org/articles/epjn/pdf/2019/01/epjn180012.pdf         #
# "geometry" takes as imput, the material created by the function             #
# "materials.py".                                                             #
#                                                                             #
# By Morten Nygaard for the DTU student project "Impact of temperature        #
# feedback on reactivity parameters in the molten salt fast reactor"          #
# under the DTU Centre of Nuclear energy technology (ongoing per June 2024).  #
###############################################################################

import openmc

def geometry(reflector, fuel_salt, b4c, fertile_salt):

    #%% DEFINE DIMENSION PARAMETERS
    
    # Please refer to the benchmark article Figure 2 for dimensions (in mm) and 
    # colors. Note that the colors differ from "geometry_plot.py". The EVOL 
    # reference MSFR is designed using cylinders, and cyliner shells in the 
    # z-direction limited by xy planes. 
    
    # Hights are defined as half-lengths, giving the z-distance from the 
    # model center (x,y,z) = (0,0,0).
    
    # Reflector (blue) outer dimensions (in cm)
    r_reflector         = 452.9/2 
    h_reflector         = 425.5/2
    
    # Fuel salt (yellow) outer dimensions 
    r_fuel_salt_out       = 225.5/2 + 50.0 + 20.0 + 23.7
    h_fuel_salt_out       = 225.5/2
    
    # b4c protection (orange) outer dimensions
    r_protection        = 225.5/2 + 50.0 + 20.0 
    h_protection        = 188.0/2
    
    # Outer dimensions of the the fertile blanket outer wall (blue)
    # This corresponds to the inner dimensions of the B4C protection (orange)
    r_blanket_wall_out    = 225.5/2 + 50.0
    h_blanket_wall_out    = 188.0/2
    
    # Outer Reflector Wall/Fertile blanket boundary 
    r_fertile_salt      = 225.5/2 + 2.0 + 46.0
    h_fertile_salt      = 184.0/2
    
    # Inner Reflector Wall/Fertile blanket boundary 
    r_blanket_wall_in    = 225.5/2 + 2.0
    h_blanket_wall_in    = 188.0/2
    
    # Inner fuel salt region (yellow) outer dimensions 
    r_fuel_salt_in       = 225.5/2
    h_fuel_salt_in       = 225.5/2
    
    #%% DEFINE BOUNDARIES
    
    # Boundaries of cylinders are difined from the inside and out, 
    # Innermost cylinder first, outermost cylinder last.
    
    # Cylinder are defined with axes parallel to the z-axis, such that the 
    # xz- an yz-planes gives he vertical cross sections, and the xy-plane gives  
    # the horizontal cross-section.
    
    # Cylinders are infinite in the z-axis, so each must be defined as half-space 
    # by using an upper and lower plane perpendicular to the z-axis
       
    # Inner fuel cylinder
    cyl_fuel_salt_in           = openmc.ZCylinder(r=r_fuel_salt_in, boundary_type='transmission')
    top_plane_fuel_salt_in     = openmc.ZPlane(z0=h_fuel_salt_in, boundary_type='transmission')
    bottom_plane_fuel_salt_in  = openmc.ZPlane(z0=-h_fuel_salt_in, boundary_type='transmission')
    
    # Inner wall/fertile salt boundary
    cyl_blanket_wall_in           = openmc.ZCylinder(r=r_blanket_wall_in, boundary_type='transmission')
    top_plane_blanket_wall_in     = openmc.ZPlane(z0=h_blanket_wall_in, boundary_type='transmission')
    bottom_plane_blanket_wall_in  = openmc.ZPlane(z0=-h_blanket_wall_in, boundary_type='transmission')
    
    # Outer wall/fertile salt boundary
    cyl_fertile_salt           = openmc.ZCylinder(r=r_fertile_salt, boundary_type='transmission')
    top_plane_fertile_salt     = openmc.ZPlane(z0=h_fertile_salt, boundary_type='transmission')
    bottom_plane_fertile_salt  = openmc.ZPlane(z0=-h_fertile_salt, boundary_type='transmission')
    
    # wall/protection boundary
    cyl_blanket_wall_out           = openmc.ZCylinder(r=r_blanket_wall_out, boundary_type='transmission')
    top_plane_blanket_wall_out     = openmc.ZPlane(z0=h_blanket_wall_out, boundary_type='transmission')
    bottom_plane_blanket_wall_out  = openmc.ZPlane(z0=-h_blanket_wall_out, boundary_type='transmission')
    
    # Protection/fuel boundary
    cyl_protection           = openmc.ZCylinder(r=r_protection, boundary_type='transmission')
    top_plane_protection     = openmc.ZPlane(z0=h_protection, boundary_type='transmission')
    bottom_plane_protection  = openmc.ZPlane(z0=-h_protection, boundary_type='transmission')
    
    # Fuel/reflector boundary
    cyl_fuel_salt_out             = openmc.ZCylinder(r=r_fuel_salt_out, boundary_type='transmission')
    top_plane_fuel_salt_out       = openmc.ZPlane(z0=h_fuel_salt_out, boundary_type='transmission')
    bottom_plane_fuel_salt_out    = openmc.ZPlane(z0=-h_fuel_salt_out, boundary_type='transmission')
    
    # Outer reflector boundary
    cyl_reflector           = openmc.ZCylinder(r=r_reflector, boundary_type='vacuum')
    top_plane_reflector     = openmc.ZPlane(z0=h_reflector, boundary_type='vacuum')
    bottom_plane_reflector  = openmc.ZPlane(z0=-h_reflector, boundary_type='vacuum')
    
    
    #%% DEFINE REGIONS
    
    # Fertile salt region (pink) CHECK
    fertile_salt_region = +cyl_blanket_wall_in & -cyl_fertile_salt \
    & -top_plane_fertile_salt & +bottom_plane_fertile_salt 
        
    # Protection region (orange) CHECK
    protection_region = +cyl_blanket_wall_out & -cyl_protection \
    & -top_plane_protection & +bottom_plane_protection 
            
    # Fuel salt region (yellow) CHECK
    fuel_salt_region = \
    (+cyl_protection & -cyl_fuel_salt_out \
    & -top_plane_fuel_salt_out & +bottom_plane_fuel_salt_out)\
    |   (-cyl_fuel_salt_out & -top_plane_fuel_salt_out \
    & +top_plane_protection)\
    |   (-cyl_fuel_salt_out & +bottom_plane_fuel_salt_out \
    & -bottom_plane_protection)\
    |   (-cyl_fuel_salt_in & -top_plane_fuel_salt_in \
    & +bottom_plane_fuel_salt_in)
    
    # Reflector and 20 mm wall region outside fertile blanket (blue) CHECK
    reflector_region = \
    (-cyl_reflector & +cyl_fuel_salt_out \
    &   -top_plane_fuel_salt_out & +bottom_plane_fuel_salt_out) \
    | (-cyl_reflector \
    & -top_plane_reflector & +top_plane_fuel_salt_out)\
    |   (-cyl_reflector \
    & +bottom_plane_reflector & -bottom_plane_fuel_salt_out)\
    |   (+cyl_fertile_salt & -cyl_blanket_wall_out \
    & -top_plane_blanket_wall_out & +bottom_plane_blanket_wall_out)\
    |   (+cyl_fuel_salt_in  & -cyl_blanket_wall_in \
    & -top_plane_blanket_wall_in & +bottom_plane_blanket_wall_in )\
    |   (+cyl_fuel_salt_in  & -cyl_blanket_wall_out \
    & +top_plane_fertile_salt & -top_plane_blanket_wall_in )\
    |   (+cyl_fuel_salt_in  & -cyl_blanket_wall_out \
    & -bottom_plane_fertile_salt & +bottom_plane_blanket_wall_in )
    
    #%% DEFINE CELLS
    
    ### OBS: Proper materials to be defined
    
    # Reflector and 20 mm wall cell, outside fertile blanket (blue)
    reflector_cell = openmc.Cell()
    reflector_cell.fill = reflector
    reflector_cell.region = reflector_region
    
    # Fuel salt cell (yellow)
    fuel_salt_cell = openmc.Cell()
    fuel_salt_cell.fill = fuel_salt
    fuel_salt_cell.region = fuel_salt_region
    
    # Protection cell (orange)
    protection_cell = openmc.Cell()
    protection_cell.fill = b4c
    protection_cell.region = protection_region
    
    # Fertile salt cell (pink)
    fertile_salt_cell = openmc.Cell()
    fertile_salt_cell.fill = fertile_salt
    fertile_salt_cell.region = fertile_salt_region
    
    #%% Define universe
    root = openmc.Universe(cells=\
    [reflector_cell,fuel_salt_cell,protection_cell,fertile_salt_cell])
    
    geometry = openmc.Geometry(root)
    geometry.export_to_xml()
        
    return geometry

################### 
### End of code ###
###################
