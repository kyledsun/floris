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
Plots streamwise velocities upstream of baseline and blockage models for varying layouts at
hub height, 3 diameters upstream of each specified wind farm layout.
"""

input_file="../../example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Set parameters for wind farm layout
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)

### Function to find nearest value in array to specified point
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# ============================================================================================
# Copy floris interface to model baseline floris model Baseline
fi2 = copy.deepcopy(fi)

# Set induciton to false for baseline model
Ind_Opts = fi2.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction']=False
fi2.IndOpts = Ind_Opts

# Creates a turbine field with m columns and n rows
m = 3; n = 1
D = fi2.floris.farm.turbines[0].rotor_diameter
layout_x = []; layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi2.reinitialize_flow_field(layout_array=[layout_x,layout_y],wind_speed=8)

# Calculate wake
fi2.calculate_wake(Ind_Opts=Ind_Opts)

# Initialize the horizontal cut at hub height
hor_plane_base = fi2.get_hor_plane(x_bounds=[-400,3500])
# Look at data at specified x location
hplane_base = hor_plane_base.df[hor_plane_base.df.x1 == find_nearest(hor_plane_base.df.x1,-3*D)]
# ============================================================================================
# Blockage

# Set induction to true to model blockage effect
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction']=True
fi.IndOpts = Ind_Opts

# Creates list of turbine layouts with m columns and n rows
m = 3
n = [1,3,5]

# Create object to store upstream velocities
inVel = []; layout = []

# Loop through and return upstream velocities for each layout
for row in n:
    D = fi.floris.farm.turbines[0].rotor_diameter
    layout_x = []
    layout_y = []
    for i in range(m):
        for j in range(row):
            layout_x.append(j*sep*D)
            layout_y.append(i*sepy*D)

    # Reinitialize flow field with new specified layout
    fi.reinitialize_flow_field(layout_array=[layout_x,layout_y],wind_speed=8)#,turbulence_intensity=0.11)

    x_grid, y_grid, z_grid = fi.floris.farm.flow_field._discretize_upstream_of_turbine_domain()#fi.floris.farm.flow_field)

    # Calculate wake
    fi.calculate_wake(Ind_Opts=Ind_Opts,points = np.array([x_grid.flatten(),y_grid.flatten(),z_grid.flatten()]))

    # Initialize the horizontal cut at hub height
    hor_plane = fi.get_hor_plane(x_bounds=[-1260,3500])
    # Look at data at specified x location
    hplane = hor_plane.df[hor_plane.df.x1 == find_nearest(hor_plane.df.x1,-3*D)]

    inVel.append(hplane)
    layout.append(str(m)+'x'+str(row)+' Blockage')

fig,ax = plt.subplots(figsize=(7,6))
ax.plot(hplane_base.x2,hplane_base.u,label='Baseline')
for i in range(len(layout)):
    print(layout[i],'minimum vel:', inVel[i].u.min())
    ax.plot(inVel[i].x2,inVel[i].u,label=layout[i])
ax.set_title('Wind Velocity 3D Upstream of Wind Farm',fontsize=22,pad=10)
ax.set_xlabel('y [m]',fontsize=16)
ax.set_ylabel('Streamwise Velocity (u)',fontsize=16)
ax.tick_params(axis='both',labelsize=14)
ax.legend(loc = 'upper right',fontsize=12)
fig.tight_layout()

# fig.savefig('../../../../../../../Documents/Blockage Effect/Paper Figures/WES Figures/UpstreamVelProfile2')

print('----------------------------------------------')
print('Time Elapsed: ', (time.time()-tstart))
plt.show()