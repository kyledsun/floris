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
Tests the difference in yaw optimized angles between varying instances of induction calculation iterations
"""

input_file="OptLayout_2x3.json"
# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file)

Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
# Sets induction to true
Ind_Opts['induction'] = True

fi.floris.farm.flow_field.Ind_Opts = Ind_Opts

print('1st Ind Opts: ', fi.floris.farm.flow_field.Ind_Opts)

sep = 5
D = fi.floris.farm.turbines[0].rotor_diameter
numTurb = 3
layout_x = []
layout_y = []
for i in range(numTurb):
    layout_y.append(0)
    layout_x.append(i*D*sep)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y])

# Set bounds for allowable wake steering
min_yaw = 0.0
max_yaw = 25.0

# Instantiate the Optimization object
yaw_opt = YawOptimization(fi, minimum_yaw_angle=min_yaw, maximum_yaw_angle=max_yaw)

# Perform optimization
yaw_angles = yaw_opt.optimize()

# Calculate wake and get horizontal plane at turbine height for yaw optimized farm field
fi.calculate_wake(yaw_angles=yaw_angles)

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_resolution=400, y_resolution=100)

# Plot and show
fig, axs = plt.subplots(nrows = 2, ncols=1)
wfct.visualization.visualize_cut_plane(hor_plane, ax=axs[0])
axs[0].set_title("1 Iteration")

base_power = fi.get_farm_power()
# =============================================================================
# Increase Iterations
# =============================================================================

# Make a copy for floris interface with induction
fi2 = copy.deepcopy(fi)

# Read in induction options from flow field class
Ind_Opts = fi2.floris.farm.flow_field.Ind_Opts

# Increase the number of iterations
Ind_Opts['nIter'] = 20

fi2.floris.farm.flow_field.Ind_Opts = Ind_Opts

print('2nd Ind_Opts: ', fi2.floris.farm.flow_field.Ind_Opts)

# Instantiate the Optimization object
yaw_ind_opt = YawOptimization(fi2, minimum_yaw_angle=min_yaw, maximum_yaw_angle=max_yaw)

# Perform optimization
yaw_ind_angles = yaw_ind_opt.optimize()

# Calculate wake and get horizontal plane at turbine height for yaw optimized farm field
fi2.calculate_wake(yaw_angles=yaw_ind_angles)

# Initialize the horizontal cut
hor_plane = fi2.get_hor_plane(x_resolution=400, y_resolution=100)

# Plot and show
wfct.visualization.visualize_cut_plane(hor_plane, ax=axs[1])
axs[1].set_title("20 Iterations")

ind_power = fi2.get_farm_power()

print("==========================================")
print('1 Iteration:', base_power)
print('20 Iterations:', ind_power)

print("==========================================")
print("Total Power Change = %.1f%%" %(100.0 * (ind_power - base_power) / base_power))

print("==========================================")
print('Turbine Powers:')
print("------------------------------------------")
print('\t\tBaseline\tIncreased Iterations')
print("------------------------------------------")
init_powers = fi.get_turbine_power()
opt_powers = fi2.get_turbine_power()
for i in range(len(init_powers)):
    print('Turbine %d: \t%.2f\t%.2f' %(i, init_powers[i],opt_powers[i]))

print("==========================================")
print('Velocities Seen By Each Turbine:')
print("------------------------------------------")
print('\t\tBaseline\tIncreased Iterations')
print("------------------------------------------")
turbine_vel = []
indTurbine_vel = []
for turbine in fi.floris.farm.turbine_map.turbines:
    turbine_vel.append(turbine.average_velocity)
for turbine in fi2.floris.farm.turbine_map.turbines:
    indTurbine_vel.append(turbine.average_velocity)

for i in range(len(turbine_vel)):
    print('Turbine %d: \t%.2f\t\t%.2f' %(i,turbine_vel[i],indTurbine_vel[i]))

print("==========================================")
print('Optimized Yaw Angles (deg): ')
print("------------------------------------------")
print('\t\tBaseline\tIncreased Iterations')
print("------------------------------------------")
for i in range(len(yaw_angles)):
    print("Turbine %d:\t%.2f\t\t%.2f " %(i,yaw_angles[i],yaw_ind_angles[i]))

plt.show()