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

# --- Plot options
bCompute=True
bColorBar=False
nStreamlines=0
U0 =6
minSpeed=0.5
maxSpeed=1.03
levelsLines=np.sort([1.05,1.0,0.99,0.98,0.95,0.9,0.5])

# ---
ny=100
nx=ny*4
resolution=Vec3(nx, ny, 2)

input_file="Layout_1x2_v2.json"

# bounds=[-4*D,14*D,-2*D-10,2*D+10,89,90] # xmin xmax .. zmin zmax
#Induction = False # Default for inclusion of the induction zone

# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts

# Calculate wake
fi.calculate_wake(Ind_Opts=Ind_Opts)

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_resolution = resolution.x1, y_resolution = resolution.x2, Ind_Opts = Ind_Opts)

# =============================================================================
# Induction
# =============================================================================
# Make a copy for floris interface with induction
fi_ind = copy.deepcopy(fi)

# Read in induction options
Ind_Opts = fi_ind.floris.farm.flow_field.Ind_Opts

# Set induction to true to model blockage effect
Ind_Opts['induction']=True

# Calculate wake
fi_ind.calculate_wake(Ind_Opts=Ind_Opts)

# Initialize the horizontal cut
hor_ind_plane = fi_ind.get_hor_plane(x_resolution = resolution.x1, y_resolution = resolution.x2, Ind_Opts = Ind_Opts)

# =============================================================================
# Plot and Show
# =============================================================================
titles=[]
titles.append('FLORIS (original)')
titles.append('FLORIS (with induction)')

fig, axs = plt.subplots(nrows=len(titles), ncols=1, sharex=True, sharey=True)
# Plot flow field without induction
wfct.visualization.visualize_cut_plane(hor_plane, ax=axs[0])
# wfct.visualization.plot_turbines_with_fi(axs[0],fi)
axs[0].set_title(titles[0])
# Plot flow field with induction
wfct.visualization.visualize_cut_plane(hor_ind_plane, ax=axs[1])
# wfct.visualization.plot_turbines_with_fi(axs[1],fi_ind)
axs[1].set_title(titles[1])

plt.show()