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

input_file="../example_induction_input.json"

# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Calculate wake
fi.calculate_wake()

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane()

hor_plane2 = copy.deepcopy(hor_plane)
hor_plane2.df.u = hor_plane2.df.u_ind

fig, ax = plt.subplots(nrows = 3, sharex=True, sharey=True)
# Plot flow field with induction
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0])
wfct.visualization.plot_turbines_with_fi(ax[0],fi)
ax[0].set_title('Flow Field with Induction')


# Plot induction
wfct.visualization.visualize_cut_plane(hor_plane2, ax=ax[1])
wfct.visualization.plot_turbines_with_fi(ax[1],fi)
ax[1].set_title('Induction Only')

ind_opts = fi.floris.farm.flow_field.Ind_Opts
ind_opts['induction'] = False
fi.floris.farm.flow_field.Ind_Opts = ind_opts

# Calculate wake
fi.calculate_wake()

# Initialize the horizontal cut
hor_plane = fi.get_hor_plane()

# Plot flow field without induction
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[2])
wfct.visualization.plot_turbines_with_fi(ax[2],fi)
ax[2].set_title('Flow Field without Induction')
fig.tight_layout()

plt.show()