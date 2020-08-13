"""
Plots velocity plots for the flow field with and without induction and induction only.
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
# contours = np.sort([1.01,0.995,0.99,0.98,0.95,0.9,0.8,0.85])
contours = np.sort([1.01,0.99,0.98,0.95,0.9,0.8,0.7,0.6,0.5])


# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)
fi.reinitialize_flow_field(layout_array=[[0],[0]])

Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['Model'] ='VC'
Ind_Opts['nIter'] = 1
Ind_Opts['Ct_test'] = True
fi.IndOpts = Ind_Opts

fi.floris.farm.flow_field.set_ct([0.95])

D=fi.floris.farm.turbines[0].rotor_diameter
bounds=[-5*D,3*D,-2*D,2*D,89,90] # xmin xmax .. zmin zmax

# Calculate wake
fi.calculate_wake(Ind_Opts=Ind_Opts)

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane(x_resolution = resolution.x1, y_resolution = resolution.x2, x_bounds = tuple(bounds[0:2]), y_bounds = tuple(bounds[2:4]))

hor_plane2 = copy.deepcopy(hor_plane)
hor_plane2.df.u = hor_plane2.df.u_ind
print(hor_plane2.df.head())

# ===============================================================================================
# Induction Only Plot
# ===============================================================================================
fig, ax = plt.subplots()
# Plot induction
wfct.visualization.visualize_cut_plane(hor_plane2, ax=ax, levels = contours)#, minSpeed = -0.25, maxSpeed = 0.25)
wfct.visualization.plot_turbines_with_fi(ax,fi)
ax.set_title('Induction Only')

# ===============================================================================================
# Induction Only, No Induction and Combined Plot
# ===============================================================================================
# fig, ax = plt.subplots(nrows = 3, sharex=True, sharey=True)
# # Plot flow field with induction
# wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0], minSpeed = minspeed, maxSpeed = maxspeed, fig=fig, cbar=True)
# wfct.visualization.plot_turbines_with_fi(ax[0],fi)
# ax[0].set_title('Flow Field with Induction')

# # Plot induction
# wfct.visualization.visualize_cut_plane(hor_plane2, ax=ax[1], minSpeed = minspeed, maxSpeed = maxspeed, fig=fig, cbar=True)
# wfct.visualization.plot_turbines_with_fi(ax[1],fi)
# ax[1].set_title('Induction Only')

# fi2 = wfct.floris_interface.FlorisInterface(input_file)
# ind_opts = fi2.floris.farm.flow_field.Ind_Opts
# ind_opts['induction'] = False
# fi2.floris.farm.flow_field.Ind_Opts = ind_opts

# # Calculate wake
# fi2.calculate_wake()

# # Initialize the horizontal cut
# hor_plane = fi2.get_hor_plane(x_resolution = resolution.x1, y_resolution = resolution.x2)

# # Plot flow field without induction
# wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[2], minSpeed = minspeed, maxSpeed = maxspeed, fig=fig, cbar=True)
# wfct.visualization.plot_turbines_with_fi(ax[2],fi2)
# ax[2].set_title('Flow Field without Induction')
# fig.tight_layout()

for turbine in fi.floris.farm.turbine_map.turbines:
    print("Turbine Ct:", turbine.Ct)

plt.show()