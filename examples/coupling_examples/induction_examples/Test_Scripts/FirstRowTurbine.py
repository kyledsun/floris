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
Returns first row turbine power as a function of number of turbine rows in single column turbine farm.
"""

input_file="../../example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
# Set induction to true to model blockage effect
Ind_Opts['induction']=True
Ind_Opts["Model"] = 'VC'
fi.IndOpts = Ind_Opts

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with m columns and n rows
m = 1
n = np.arange(1,26)

D = fi.floris.farm.turbines[0].rotor_diameter

firstrowturb = []
for numrow in n:

    layout_x = []
    layout_y = []
    for i in range(m):
        for j in range(numrow):
            layout_x.append(j*sep*D)
            layout_y.append(i*sepy*D)

    # Reinitialize flow field with new specified layout
    fi.reinitialize_flow_field(layout_array=[layout_x,layout_y],wind_speed=8)

    x_grid, y_grid, z_grid = fi.floris.farm.flow_field._discretize_upstream_of_turbine_domain()#fi.floris.farm.flow_field)
    # Calculate wake
    fi.calculate_wake(Ind_Opts=Ind_Opts,points = np.array([x_grid.flatten(),y_grid.flatten(),z_grid.flatten()]))

    # Return Turbine Powers
    turbpow = np.array(fi.get_turbine_power())/1000    
    firstrowturb.append(turbpow[0])

print(firstrowturb)

fig,ax = plt.subplots(figsize=(7,6))
ax.plot(n,firstrowturb,marker='o',markersize=4)
ax.set_xlabel('Number of Turbine Rows',fontsize=16)
ax.set_ylabel('Power (kW)',fontsize=16)
ax.set_title('First Row Turbine Power',fontsize=22,pad=10)
ax.tick_params(axis='both',labelsize=14)
fig.tight_layout()

# fig.savefig('../../../../../../../Documents/Blockage Effect/Paper Figures/WES Figures/FirstRowTurbinePow2')

print('----------------------------------------------')
print('Time Elapsed: ', (time.time()-tstart))
plt.show()