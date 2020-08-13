import os
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
# --- Local
from floris.tools.optimization.scipy.yaw import YawOptimization
import floris.tools as wfct
from floris.utilities import Vec3
"""
Compares the turbine parameters of a wind farm with and without modelling blockage effects in the 
induction zone of turbines.
"""

input_file="../OptLayout_2x3.json"
# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file)

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with n rows and m columns
n = 5
m = 1

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y])

# Calculate wake and get horizontal plane at turbine height for yaw optimized farm field
fi.calculate_wake()

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_resolution=400, y_resolution=100, y_bounds=[-500,500])

# Plot and show
fig, axs = plt.subplots(nrows = 2, ncols=1)
wfct.visualization.visualize_cut_plane(hor_plane, ax=axs[0])
axs[0].set_title("No Induction")

base_power = fi.get_farm_power()
# =============================================================================
# Induction
# =============================================================================

# Make a copy for floris interface with induction
fi_ind = copy.deepcopy(fi)

# Read in induction options from flow field class
Ind_Opts = fi_ind.floris.farm.flow_field.Ind_Opts

# Sets induction to true
Ind_Opts['induction'] = True
fi_ind.floris.farm.flow_field.Ind_Opts = Ind_Opts

# Calculate wake and get horizontal plane at turbine height for yaw optimized farm field
fi_ind.calculate_wake()

# Initialize the horizontal cut
hor_plane = fi_ind.get_hor_plane(x_resolution=400, y_resolution=100, y_bounds=[-500,500])

# Plot and show
wfct.visualization.visualize_cut_plane(hor_plane, ax=axs[1])
axs[1].set_title("With Blockage Effect")

ind_power = fi_ind.get_farm_power()

print("================================================================")
print('Layout: %dx%d, Induction Iterations: %d, Streamwise Separation: %d' %(m,n,Ind_Opts['nIter'],sep))
print("================================================================")
print("Farm Powers: ")
print("----------------------------------------------------------------")
print('No Induction Power:', base_power)
print('Induction Power:', ind_power)

print("================================================================")
print("Total Power Change = %.2f%%" %(100.0 * (ind_power - base_power) / base_power))

print("================================================================")
print('Turbine Powers:')
print("----------------------------------------------------------------")
print('\t\tNo Induction\tInduction\tDifference')
print("----------------------------------------------------------------")
init_powers = fi.get_turbine_power()
opt_powers = fi_ind.get_turbine_power()
for i in range(len(init_powers)):
    print('Turbine %d: \t%.2f\t%.2f\t%.2f' %(i, init_powers[i],opt_powers[i],(opt_powers[i] - init_powers[i])))

print("================================================================")
print('Velocities Seen By Each Turbine:')
print("----------------------------------------------------------------")
print('\t\tNo Induction\tInduction\tDifference')
print("----------------------------------------------------------------")
turbine_vel = []
indTurbine_vel = []
for turbine in fi.floris.farm.turbine_map.turbines:
    turbine_vel.append(turbine.average_velocity)
for turbine in fi_ind.floris.farm.turbine_map.turbines:
    indTurbine_vel.append(turbine.average_velocity)

for i in range(len(turbine_vel)):
    print('Turbine %d: \t%.4f\t\t%.4f\t\t%.4f' %(i,turbine_vel[i],indTurbine_vel[i],(indTurbine_vel[i]-turbine_vel[i])))

plt.show()