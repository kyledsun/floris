import numpy as np
import copy
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
# --- Local
from floris.tools.optimization.scipy.yaw import YawOptimization
import floris.tools as wfct
from floris.utilities import Vec3

"""
Compares the power outputs and velocities seen by each turbine in a 2x3 wind farm for
non-optimized and yaw-optimized instances.
- Can be set to include blockage effects
"""

# --- Resolution Parameters
ny= 100
nx=ny*4
resolution=Vec3(nx, ny, 2)

input_file="OptLayout_2x3.json"
# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file)

D = fi.floris.farm.turbines[0].rotor_diameter
bounds=[-4*D,16*D,-2*D-10,5*D+10,89,90] # xmin xmax .. zmin zmax

# Make a copy for optimization instance
fi_opt = copy.deepcopy(fi)

# Read in induction options from input file
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts

Ind_Opts['induction'] = True

fi.IndOpts = Ind_Opts
fi_opt.IndOpts = Ind_Opts

# Calculate wake and get horizontal plane at turbine height for original farm field
fi.calculate_wake()
hor_plane = fi.get_hor_plane(
        x_resolution = resolution.x1,
        y_resolution = resolution.x2,
        # x_bounds = tuple(bounds[0:2]),
        # y_bounds = tuple(bounds[2:4]),
        Ind_Opts = Ind_Opts
        )

# Store value for total wind farm power
power_initial = fi.get_farm_power()

# fig, ax = plt.subplots()
# wfct.visualization.visualize_cut_plane(hor_plane, ax=ax)
# # wfct.visualization.plot_turbines_with_fi(ax,fi_opt)
# ax.set_title("Baseline Case for U = 8 m/s, Wind Direction = 270$^\circ$")

# ==================================================================================
# Yaw angle optimization
# ==================================================================================
# Set bounds for allowable wake steering
min_yaw = 0.0
max_yaw = 25.0

# Instantiate the Optimization object
yaw_opt = YawOptimization(fi_opt, minimum_yaw_angle=min_yaw, maximum_yaw_angle=max_yaw)

# Perform Optimization
yaw_angles = yaw_opt.optimize()

# Calculate wake and get horizontal plane at turbine height for yaw optimized farm field
fi_opt.calculate_wake(Ind_Opts=Ind_Opts,yaw_angles=yaw_angles)

hor_plane_opt = fi_opt.get_hor_plane(
        x_resolution = resolution.x1,
        y_resolution = resolution.x2,
        # x_bounds = tuple(bounds[0:2]),
        # y_bounds = tuple(bounds[2:4]),
        Ind_Opts = Ind_Opts
        )

power_opt = fi_opt.get_farm_power()

# ==================================================================================
# Plot Results
# ==================================================================================

# --- Plot and show
fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(12.0, 6.0))
wfct.visualization.visualize_cut_plane(hor_plane, ax=axs[0])
wfct.visualization.plot_turbines_with_fi(axs[0],fi)
axs[0].set_title("Baseline Case for U = 8 m/s, Wind Direction = 270$^\circ$")

wfct.visualization.visualize_cut_plane(hor_plane_opt, ax=axs[1])
wfct.visualization.plot_turbines_with_fi(axs[1],fi_opt)
axs[1].set_title("Optimized Case for U = 8 m/s, Wind Direction = 270$^\circ$")

print("==========================================")
print("optimized yaw angles = ")
for i in range(len(yaw_angles)):
    print("Turbine ", i, "=", yaw_angles[i], " deg")

print("==========================================")
print('Original Power:', power_initial)
print('Optimized Power:', power_opt)

print("==========================================")
print("Total Power Gain = %.1f%%" %(100.0 * (power_opt - power_initial) / power_initial))

print("==========================================")
print('Turbine Powers:')
print("------------------------------------------")
print('\t\tBaseline\tOptimized')
print("------------------------------------------")
init_powers = fi.get_turbine_power()
opt_powers = fi_opt.get_turbine_power()
for i in range(len(init_powers)):
    print('Turbine %d: \t%.2f\t%.2f' %(i, init_powers[i],opt_powers[i]))

print("==========================================")
print('Velocities Seen By Each Turbine:')
print("------------------------------------------")
print('\t\tBaseline\tOptimized')
print("------------------------------------------")
turbine_vel = []
optTurbine_vel = []
for turbine in fi.floris.farm.turbine_map.turbines:
    turbine_vel.append(turbine.average_velocity)
for turbine in fi_opt.floris.farm.turbine_map.turbines:
    optTurbine_vel.append(turbine.average_velocity)

for i in range(len(turbine_vel)):
    print('Turbine %d: \t%.2f\t\t%.2f' %(i,turbine_vel[i],optTurbine_vel[i]))

plt.show()