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
Plots horizontal plane and upstream hub height wind velocities at x = -252. Returns farm and turbine powers.
- Can also set blockcomp flag to True to compare plots with and without blockages.
"""

blockcomp = False

input_file="../../../example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Read in upstream velocity profile from SOWFA case.
InputSpeeds = pd.read_csv('../../../../sowfa_comparisons/SowfaVelocities/SowfaInputVel3x1.csv')

def inSpeeds(numpoints,df):
    locations = np.array([float(i) for i in df[df.NumPoints == numpoints].Locations.values[0].strip('][').split(', ')]) - 320
    velocities = np.array([float(i) for i in df[df.NumPoints == numpoints].Velocities.values[0].strip('][').split(', ')]) #* (8.38 / 8)
    return locations,velocities

# Return wind locations and velocities from sowfa velocity files at sampled points
samploc,sampvel = inSpeeds(3,InputSpeeds)
xwind = np.ones(len(samploc))-252

levelsLines=np.arange(0.97,1.01,0.005) * fi.floris.farm.wind_map.input_speed[0]

# print('Locations: ', samploc)
# print('Velocities: ', sampvel)

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

# print([xwind,list(samploc)])

# Reinitialize flow field with new specified layout
# fi.reinitialize_flow_field(layout_array=[layout_x,layout_y])
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y], wind_speed=list(sampvel+[0.0559,0.5892,0.2019]),wind_layout=[xwind,list(samploc)],turbulence_intensity=0.11)

# print('Wind Vel: ', sampvel+ [0.0559,0.5892,0.2019])

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
# Set induction to true to model blockage effect
Ind_Opts['induction']=False
Ind_Opts["Model"] = 'VC'
fi.IndOpts = Ind_Opts

# Calculate wake
fi.calculate_wake()

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_resolution=400, y_resolution=100)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

hplane = hor_plane.df[hor_plane.df.x1 == find_nearest(hor_plane.df.x1,-252)]

# Plot and Show Horizontal Plane
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax, fig=fig, cbar=True)#,minSpeed=4,maxSpeed=9)
wfct.visualization.plot_turbines_with_fi(ax,fi)
ax.set_title('FLORIS Horizontal Plane',fontsize=20)

if blockcomp:

    fi2 = copy.deepcopy(fi)
    Ind_Opts = fi2.floris.farm.flow_field.Ind_Opts
    Ind_Opts['induction']=True
    fi2.IndOpts = Ind_Opts

    fi2.calculate_wake()

    # Initialize the horizontal cut
    hor_plane2 = fi2.get_hor_plane(x_resolution=400, y_resolution=100)
    hplane2 = hor_plane2.df[hor_plane2.df.x1 == find_nearest(hor_plane2.df.x1,-252)]
    # # Initialize the cross plane cut at x = -252 (SOWFA input conditions)
    # cross_plane2 = fi2.get_cross_plane(-252,y_resolution=400)
    # cplane2 = cross_plane2.df[cross_plane2.df.x2 == find_nearest(cross_plane2.df.x2,90.0)]

    # Plot and Show Horizontal Plane
    fig, ax = plt.subplots()
    wfct.visualization.visualize_cut_plane(hor_plane2, ax=ax, fig=fig, cbar=True)#,minSpeed=4,maxSpeed=9)
    wfct.visualization.plot_turbines_with_fi(ax,fi2)
    ax.set_title('FLORIS Horizontal Plane Block',fontsize=20)

    # Plot Comparing horizontal velocity planes
    fig,ax = plt.subplots(1,2,figsize=(12,5))
    wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0], fig=fig, cbar=True)#,minSpeed=4,maxSpeed=9)
    wfct.visualization.plot_turbines_with_fi(ax[0],fi)
    ax[0].set_title('FLORIS Horizontal Plane',fontsize=22)
    ax[0].set_xlabel('x [m]',fontsize=16); ax[0].set_ylabel('y [m]',fontsize=16) 
    wfct.visualization.visualize_cut_plane(hor_plane2, ax=ax[1], fig=fig, cbar=True)#,minSpeed=4,maxSpeed=9)
    wfct.visualization.plot_turbines_with_fi(ax[1],fi2)
    ax[1].set_title('FLORIS Horizontal Plane Block',fontsize=22)
    ax[1].set_xlabel('x [m]',fontsize=16); ax[1].set_ylabel('y [m]',fontsize=16) 
    fig.tight_layout()

    # Plot streamwise velocities at hub height at input
    fig,ax = plt.subplots()
    ax.plot(hplane.x2, hplane.u, label = 'baseline')
    ax.plot(hplane2.x2, hplane2.u, label = 'blockage')
    ax.set_xlabel('y [m]', fontsize=14)
    ax.set_ylabel('stream wise velocity (u)',fontsize=14)
    ax.set_title('Wind Speed at Hub Height at x=-252',fontsize=18)
    ax.legend(fontsize=12)

# Plot streamwise velocities at hub height at input
fig,ax = plt.subplots()
ax.plot(hplane.x2, hplane.u)
# ax.scatter(hplane.x1, hplane.u, s= 5)
ax.set_xlabel('y [m]',fontsize=16)
ax.set_ylabel('Streamwise Velocity (u)',fontsize=16)
ax.set_title('Heterogeneous Input Velocity',fontsize=24)
ax.tick_params(axis='both',labelsize=12)
fig.tight_layout()

print('==============================================')
print("Farm Power: ", fi.get_farm_power())
print('==============================================')
print("Turbine Powers:")
print('----------------------------------------------')
turbpow = np.array(fi.get_turbine_power())/1000
turbine_powers = pd.DataFrame()
for i in range(m):
    temp = []
    for j in range(n):
        temp.append(turbpow[i*n+j])
    turbine_powers['Column_'+str(i)]=temp
print(turbine_powers)
print('----------------------------------------------')
print('Time Elapsed: ', (time.time()-tstart))

plt.show()