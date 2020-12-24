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
Plots horizontal plane with blockage effects. Returns farm and turbine powers.
"""

input_file="example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with m columns and n rows
m = 1
n = 2

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y], wind_speed=8)

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
# Set induction to true to model blockage effect
Ind_Opts['induction']=True
fi.IndOpts = Ind_Opts

# Calculate wake
fi.calculate_wake()

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane()

print('==============================================')
print("Farm Power (kW): ", fi.get_farm_power()/1000)
print('==============================================')
print("Turbine Powers (kW):")
print('----------------------------------------------')
turbpow = np.array(fi.get_turbine_power())/1000 # Returns turbine powers and converts to kW
turbine_powers = pd.DataFrame()
for i in range(m):
    temp = []
    for j in range(n):
        temp.append(turbpow[i*n+j])
    turbine_powers['Column_'+str(i)]=temp
print(turbine_powers)

print('==============================================')
if fi.IndOpts['induction']:
    levelsLines=np.arange(0.97,1.01,0.005) * fi.floris.farm.wind_map.input_speed[0]
else: levelsLines = None

# Plot and Show
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax, fig=fig, cbar=True, levels=levelsLines)
wfct.visualization.plot_turbines_with_fi(ax,fi)
if Ind_Opts['induction']:
    ax.set_title('Floris with Blockage',fontsize=14)
else:
    ax.set_title('Baseline Floris',fontsize=14)
fig.tight_layout()

print('Time Elapsed: ', (time.time()-tstart))
plt.show()