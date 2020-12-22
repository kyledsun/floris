import numpy as np
import pandas as pd
import copy
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
# --- Local
import floris.tools as wfct
from floris.utilities import Vec3
import time
tstart = time.time()
"""
Plots horizontal plane and upstream hub height wind velocities and returns farm and turbine powers for baseline and blockage FLORIS
"""

input_file="../../../example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

InputSpeeds = pd.read_csv('../../../../sowfa_comparisons/SowfaVelocities/SowfaInputVel3x1.csv')

def inSpeeds(numpoints,df):
    locations = np.array([float(i) for i in df[df.NumPoints == numpoints].Locations.values[0].strip('][').split(', ')]) - 320
    velocities = np.array([float(i) for i in df[df.NumPoints == numpoints].Velocities.values[0].strip('][').split(', ')]) #* (8.38 / 8)
    return locations,velocities

# Return wind locations and velocities from sowfa velocity files at sampled points
samploc,sampvel = inSpeeds(3,InputSpeeds)
xwind = np.ones(len(samploc))-252

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
# Set induction to true to model blockage effect
Ind_Opts['induction']=True
fi.IndOpts = Ind_Opts

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with m columns and n rows
m = 3
n = 1

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y],wind_speed=8.18,turbulence_intensity=0.11)

if fi.floris.farm.flow_field.Ind_Opts['induction']:
    levelsLines=np.arange(0.97,1.01,0.005) * fi.floris.farm.wind_map.input_speed[0]
else: levelsLines = None
levelsLines = None

# Reinitialize flow field with heterogeneous velocity input
fi_het = copy.deepcopy(fi)
fi_het.reinitialize_flow_field(wind_speed=list(sampvel+[0.0559,0.5892,0.2019]),wind_layout=[xwind,list(samploc)],turbulence_intensity=0.11)

# Calculate wake
fi.calculate_wake()
fi_het.calculate_wake()

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_resolution=400, y_resolution=100, x_bounds=[-250,4250])
hor_plane_het = fi_het.get_hor_plane(x_resolution=400, y_resolution=100,x_bounds=[-250,4250])

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

hplane = hor_plane.df[hor_plane.df.x1 == find_nearest(hor_plane.df.x1,-252)]
hplane_het = hor_plane_het.df[hor_plane_het.df.x1 == find_nearest(hor_plane_het.df.x1,-252)]

# Plot and Show Horizontal Plane
# fig,ax = plt.subplots(1,2,figsize=(12,5))
fig,ax = plt.subplots(2,1,figsize=(7,6))
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0], fig=fig, cbar=True, levels=levelsLines,minSpeed=4,maxSpeed=9)
wfct.visualization.plot_turbines_with_fi(ax[0],fi)
ax[0].set_title('FLORIS Homogeneous Input Velocity',fontsize=20,pad=15)
# ax[0].set_xlabel('x [m]',fontsize=16); ax[0].set_ylabel('y [m]',fontsize=16)
ax[0].tick_params(axis='both',labelsize=12)

wfct.visualization.visualize_cut_plane(hor_plane_het, ax=ax[1], fig=fig, cbar=True, levels=levelsLines,minSpeed=4,maxSpeed=9)
wfct.visualization.plot_turbines_with_fi(ax[1],fi_het)
ax[1].set_title('FLORIS Heterogeneous Input Velocity',fontsize=20,pad=15)
# ax[1].set_xlabel('x [m]',fontsize=16); ax[1].set_ylabel('y [m]',fontsize=16)
ax[1].tick_params(axis='both',labelsize=12)
fig.tight_layout()

# Plot streamwise velocities at hub height at input
fig,ax = plt.subplots()
ax.plot(hplane.x2, hplane.u, label = 'Homogeneous')
ax.plot(hplane_het.x2, hplane_het.u, label = 'Heterogeneous')
ax.set_xlabel('y [m]',fontsize=14)
ax.set_ylabel('Streamwise Velocity (u)',fontsize=16)
ax.set_title('Heterogeneous Input Velocity',fontsize=24,pad=15)
ax.tick_params(axis='both',labelsize=12)
ax.legend(fontsize=12)
fig.tight_layout()

print('==============================================')
print("Homogeneous Farm Power: ", fi.get_farm_power())
print("Heterogeneous Farm Power: ", fi_het.get_farm_power())
print('----------------------------------------------')
print('Time Elapsed: ', (time.time()-tstart))

plt.show()