
############################################################################
# Generates cross-sectional plots of the model in the available .xml files #
############################################################################

import openmc

#plotting geometry
plots = openmc.Plots()

x_width = 500
y_width = 500
res = 1000

#xy plot
p1 = openmc.Plot()
#p1.origin=[1,1,1]
p1.basis = 'xy'
p1.width = (x_width,y_width)
p1.pixels = (res,res)
p1.filename='plot_xy'
p1.color_by = 'material'

#xz plot
p2 = openmc.Plot()
#p2.origin=[1,1,1]
p2.basis = 'xz'
p2.width = (x_width,y_width)
p2.pixels = (res,res)
p2.filename='plot_xz'
p2.color_by = 'material'

p3 = openmc.Plot()
#p3.origin=[1,1,1]
p3.basis = 'yz'
p3.width = (x_width,y_width)
p3.pixels = (res,res)
p3.filename='plot_yz'
p3.color_by = 'material'

plots.append(p1)
plots.append(p2)
plots.append(p3)
plots.export_to_xml()

openmc.plot_geometry()

###################
### End of code ###
###################
