"""
Plots induction velocity plots with and without ground effect.
"""

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

input_file="../example_induction_input.json"
# --- Resolution Parameters
ny= 100
nx=ny*4
resolution=Vec3(nx, ny, 2)

minspeed = None #2
maxspeed = None #7

# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)
fi.reinitialize_flow_field(layout_array=[[0],[0]])
fi2 = copy.deepcopy(fi)

Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['Model'] ='VC'
Ind_Opts['nIter'] = 1
Ind_Opts['Ground'] = True
fi.IndOpts = Ind_Opts

# Calculate wake
fi.calculate_wake(Ind_Opts=Ind_Opts)
# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_resolution = resolution.x1, y_resolution = resolution.x2)
print(hor_plane.df.head())

hor_plane2 = copy.deepcopy(hor_plane)
hor_plane2.df.u = hor_plane2.df.u_ind

horplane = fi.get_hor_plane(-fi.floris.farm.turbines[0].hub_height, x_resolution = resolution.x1, y_resolution = resolution.x2)
horplane2 = copy.deepcopy(horplane)
horplane2.df.u = horplane2.df.u_ind

# Vertical Planes
y_plane = fi.get_y_plane(y_loc=0,x_resolution = resolution.x1, z_resolution = resolution.x2)
y_plane2 = copy.deepcopy(y_plane)
y_plane2.df.u = y_plane2.df.u_ind
y_plane3 = copy.deepcopy(y_plane)
y_plane3.df.u = y_plane3.df.w_ind

# Horizontal Plots
fig, ax = plt.subplots(nrows=2,sharex=True,sharey=True)
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0])
ax[0].set_title('Horizontal Plot at Hub Height')
wfct.visualization.visualize_cut_plane(horplane, ax=ax[1])
ax[1].set_title('Horizontal Plot at (-) Hub Height')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))

# Streamwise Induction only Horizontal Plots
fig, ax = plt.subplots(nrows=2,sharex=True,sharey=True)
wfct.visualization.visualize_cut_plane(hor_plane2, ax=ax[0])
ax[0].set_title('Horizontal Plot at Hub Height (u_ind)')
wfct.visualization.visualize_cut_plane(horplane2, ax=ax[1])
ax[1].set_title('Horizontal Plot at (-) Hub Height (u_ind)')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))

# Vertical Plane Plots
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(y_plane, ax=ax)
ax.set_title('Vertical Plot at Y=0')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))

# Induction only Vertical Plots
fig, ax = plt.subplots(nrows=2,sharex=True,sharey=True)
wfct.visualization.visualize_cut_plane(y_plane2, ax=ax[0])
ax[0].set_title('Vertical Plot at y = 0 (u_ind)')
wfct.visualization.visualize_cut_plane(y_plane3, ax=ax[1])
ax[1].set_title('Vertical Plot at y = 0 (w_ind)')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))

# ====================================================================
# Ground Effect Removed
# ====================================================================
Ind_Opts['Ground'] = False
fi2.IndOpts = Ind_Opts

# Calculate wake
fi2.calculate_wake(Ind_Opts=Ind_Opts)
# Initialize the horizontal cut
hor_plane = fi2.get_hor_plane(x_resolution = resolution.x1, y_resolution = resolution.x2)
print(hor_plane.df.head())

hor_plane2 = copy.deepcopy(hor_plane)
hor_plane2.df.u = hor_plane2.df.u_ind

horplane = fi2.get_hor_plane(-fi2.floris.farm.turbines[0].hub_height, x_resolution = resolution.x1, y_resolution = resolution.x2)
horplane2 = copy.deepcopy(horplane)
horplane2.df.u = horplane2.df.u_ind

# Vertical Planes
y_plane = fi2.get_y_plane(y_loc=0,x_resolution = resolution.x1, z_resolution = resolution.x2)
y_plane2 = copy.deepcopy(y_plane)
y_plane2.df.u = y_plane2.df.u_ind
y_plane3 = copy.deepcopy(y_plane)
y_plane3.df.u = y_plane3.df.w_ind

# Horizontal Plots
fig, ax = plt.subplots(nrows=2,sharex=True,sharey=True)
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0])
ax[0].set_title('Horizontal Plot at Hub Height')
wfct.visualization.visualize_cut_plane(horplane, ax=ax[1])
ax[1].set_title('Horizontal Plot at (-) Hub Height')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))

# Streamwise Induction only Horizontal Plots
fig, ax = plt.subplots(nrows=2,sharex=True,sharey=True)
wfct.visualization.visualize_cut_plane(hor_plane2, ax=ax[0])
ax[0].set_title('Horizontal Plot at Hub Height (u_ind)')
wfct.visualization.visualize_cut_plane(horplane2, ax=ax[1])
ax[1].set_title('Horizontal Plot at (-) Hub Height (u_ind)')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))

# Vertical Plane Plots
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(y_plane, ax=ax)
ax.set_title('Vertical Plot at Y=0')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))

# Induction only Vertical Plots
fig, ax = plt.subplots(nrows=2,sharex=True,sharey=True)
wfct.visualization.visualize_cut_plane(y_plane2, ax=ax[0])
ax[0].set_title('Vertical Plot at y = 0 (u_ind)')
wfct.visualization.visualize_cut_plane(y_plane3, ax=ax[1])
ax[1].set_title('Vertical Plot at y = 0 (w_ind)')
fig.suptitle('Ground Effect: '+str(Ind_Opts['Ground']))
fig.tight_layout()

plt.show()