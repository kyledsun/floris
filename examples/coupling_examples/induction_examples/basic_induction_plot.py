import numpy as np
import copy
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
# --- Local
import floris.tools as wfct
from floris.utilities import Vec3

"""
Plots horizontal plane with blockage effects
"""

ny=100
nx=ny*4
resolution=Vec3(nx, ny, 2)

input_file="../Layout_1x2.json"

# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

D = fi.floris.farm.turbines[0].rotor_diameter
bounds=[-4*D,16*D,-2*D-10,2*D+10,89,90] # xmin xmax .. zmin zmax

fi.reinitialize_flow_field(layout_array=[[0],[0]])

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
# Set induction to true to model blockage effect
Ind_Opts['induction']=True
Ind_Opts["Model"] = 'VC'

fi.IndOpts = Ind_Opts

# Calculate wake
fi.calculate_wake(yaw_angles=[20], Ind_Opts=Ind_Opts)

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(
    # x_resolution = resolution.x1,
    # y_resolution = resolution.x2,
    x_bounds = tuple(bounds[0:2]),
    y_bounds = tuple(bounds[2:4]),
    Ind_Opts = Ind_Opts)

# Plot and Show
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax)
wfct.visualization.plot_turbines_with_fi(ax,fi)
ax.set_title('Floris with Induction')
plt.show()