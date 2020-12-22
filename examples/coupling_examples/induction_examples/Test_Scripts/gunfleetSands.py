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
Plots horizontal plane with blockage effects and normalized first row turbine powers for Gunfleet Sands Wind Farm.
"""

input_file="../../example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Set paramters for iteration test
sep = 3.5 # streamwise separation for turbines (*D)
sepy = 7.5 # spanwise spearation between turbines (*D)
# sep = 2.8 # streamwise separation for turbines (*D)
# sepy = 5.75 # spanwise spearation between turbines (*D)
# Creates a turbine field with m columns and n rows
m = 6
n = 7

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)
        
layoutx = [sep*D*0,sep*D*0,sep*D*0,sep*D*1,sep*D*1,sep*D*1]
layouty = [-sepy*D*1,-sepy*D*2,-sepy*D*3,-sepy*D*1,-sepy*D*2,-sepy*D*3]

layout_x = layout_x + layoutx
layout_y = layout_y + layouty

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y],wind_speed=8.0,turbulence_intensity=0.06,wind_direction=235)

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
# Set induction to true to model blockage effect
Ind_Opts['induction']=True
Ind_Opts["Model"] = 'VC'
fi.IndOpts = Ind_Opts

# Calculate wake
fi.calculate_wake(Ind_Opts=Ind_Opts)

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_bounds=[-1500,4000],y_bounds=[-3500,5500],y_resolution=400)

print('==============================================')
print("Farm Power (kW): ", fi.get_farm_power()/1000)
print('==============================================')
print("Turbine Powers:")
turbpow = np.array(fi.get_turbine_power())/1000 # Returns turbine powers and converts to kW
print(turbpow)
print('==============================================')
df = pd.DataFrame(columns=['x','y','powers'])
df['x'] = layout_x
df['y'] = layout_y
df['powers'] = turbpow

firstRow = df[df.x==0].sort_values(['y'],ascending=False).reset_index().drop(columns = ['index'])
firstRow.index = firstRow.index+1
print(firstRow)

fig,ax = plt.subplots()
ax.plot(firstRow.index,firstRow.powers/firstRow.powers[1],marker='o')
ax.set_xlabel('Turbine Number')
ax.set_ylabel('Pn/P1')
ax.set_title('Gunfleet Sands Normalized Power',fontsize=14)

fig,ax = plt.subplots(figsize=(8,4))
ax.plot(firstRow.index,firstRow.powers/firstRow.powers[1],marker='o')
ax.set_xlabel('Turbine Number')
ax.set_ylabel('Pn/P1')
ax.set_title('Gunfleet Sands Normalized Power',fontsize=14)
ax.set_yticks(np.arange(0.9,1.1,0.05))
ax.grid('both')

if fi.floris.farm.flow_field.Ind_Opts['induction']:
    levelsLines=np.arange(0.97,1.01,0.005) * fi.floris.farm.wind_map.input_speed[0]
else: levelsLines = None

# Plot and Show
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax, fig=fig, cbar=True,levels=levelsLines)#,minSpeed=0,maxSpeed=8.0)
wfct.visualization.plot_turbines_with_fi(ax,fi)
ax.set_title('WS: %.2f, TI: %.2f, WD: %d' %(fi.floris.farm.wind_map.input_speed[0],fi.floris.farm.wind_map.input_ti[0],fi.floris.farm.wind_map.input_direction[0]))
# if Ind_Opts['induction']:
#     ax.set_title('Floris with Blockage',fontsize=14)
# else:
#     ax.set_title('Baseline Floris',fontsize=14)

print('==============================================')
print('Time Elapsed: ', (time.time()-tstart))
plt.show()