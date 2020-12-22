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
Compares wind farm powers with and without modelling blockage for various layouts (Poster Figure).
"""

Horizontal_Plots = True
TurbineComp = True
input_file="../example_induction_input.json"
# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file)

# Set paramters for iteration test
D = fi.floris.farm.turbines[0].rotor_diameter
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)

# Creates a turbine field with m columns and n rows
m = [1,1,3,3,5]
n = [3,5,3,5,5]

layout = []
base_farm_pow = []
block_farm_pow = []
for i in range(len(n)):
    layout.append(str(m[i])+'x'+str(n[i]))

    layout_x = []
    layout_y = []
    for j in range(m[i]):
        for k in range(n[i]):
            layout_x.append(k*sep*D)
            layout_y.append(j*sepy*D)

    # Reinitialize flow field with new specified layout
    fi.reinitialize_flow_field(layout_array=[layout_x,layout_y],wind_speed=8)
    
    # Set induction to false to model each layout using baseline FLORIS
    Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
    Ind_Opts['induction'] = False
    fi.floris.farm.flow_field.Ind_Opts = Ind_Opts

    # Calculate wake and get horizontal plane at turbine height for yaw optimized farm field
    fi.calculate_wake()

    base_farm_pow.append(fi.get_farm_power()/1000)
    # ---------------------------------------------------------------------------------
    # Make a copy for floris interface with induction
    fi_ind = copy.deepcopy(fi)

    # Read in induction options from flow field class
    Ind_Opts = fi_ind.floris.farm.flow_field.Ind_Opts
    Ind_Opts['induction'] = True
    fi_ind.floris.farm.flow_field.Ind_Opts = Ind_Opts

    # Calculate wake and get horizontal plane at turbine height for yaw optimized farm field
    fi_ind.calculate_wake()

    if Horizontal_Plots:
        # Initialize the horizontal cut
        hor_plane = fi.get_hor_plane(x_resolution=400, y_resolution=100)

        # Plot and show
        fig, axs = plt.subplots(nrows = 2, ncols=1)
        wfct.visualization.visualize_cut_plane(hor_plane, ax=axs[0], fig=fig, cbar=True)
        wfct.visualization.plot_turbines_with_fi(axs[0],fi)
        axs[0].set_title("FLORIS Baseline", fontsize=16)

        # Initialize the horizontal cut
        hor_plane2 = fi_ind.get_hor_plane(x_resolution=400, y_resolution=100)

        # Plot and show
        wfct.visualization.visualize_cut_plane(hor_plane2, ax=axs[1], fig=fig, cbar=True)
        wfct.visualization.plot_turbines_with_fi(axs[1],fi)
        axs[1].set_title("FLORIS with Blockage Effects", fontsize=16)
        fig.tight_layout()

    if TurbineComp:
        BasePower = np.array(fi.get_turbine_power())
        BlockPower = np.array(fi_ind.get_turbine_power())

        # Plot comparing individual turbine powers
        fig,ax = plt.subplots()
        ax.plot(np.arange(len(BasePower)),(BlockPower-BasePower)/BasePower*100)
        ax.plot(np.arange(len(BasePower)),np.zeros(len(BasePower)),linestyle='dashed',linewidth=0.75,color = 'black')
        ax.scatter(np.arange(len(BasePower)),(BlockPower-BasePower)/BasePower*100)
        ax.tick_params(axis='y',labelsize=12)
        ax.set_xticks(np.arange(len(BlockPower)))
        ax.set_xticklabels(['T'+str(i) for i in range(len(BlockPower))],fontsize=14)
        ax.set_ylabel('Relative Power Change \nDue to Blockage (%)',fontsize=14)
        ax.set_xlabel('Turbine',fontsize=14)
        ax.set_title('%dx%d Turbine Blockage Comparison' %(m[i],n[i]) ,fontsize=20)
        fig.tight_layout()

    block_farm_pow.append(fi_ind.get_farm_power()/1000)

base_farm_pow = np.array(base_farm_pow)
block_farm_pow = np.array(block_farm_pow)
Powerdiff = base_farm_pow-block_farm_pow
Percerr = (block_farm_pow - base_farm_pow)/base_farm_pow*100
Percdiff = abs(block_farm_pow - base_farm_pow)/((block_farm_pow + base_farm_pow) / 2) * 100

# Bar plot comparing baseline and blockage farm powers
x = np.arange(len(block_farm_pow))*0.55
width = 0.2

fig,ax = plt.subplots()
bar1 = ax.bar(x-width/2,base_farm_pow,width,label='Baseline')
bar2 = ax.bar(x+width/2,block_farm_pow,width,label='FLORIS w/ Blockage',color='goldenrod')
ax.set_xticks(x)
ax.set_xticklabels(layout,fontsize=14)
ax.set_ylabel('Farm Powers (kW)',fontsize=14)
ax.set_xlabel('Layout',fontsize=14)
ax.set_ylim([0,20000])

# Add counts above the two bar graphs
for i,rect in enumerate(bar1):
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width(), height+250, '%.2f%%' % Percerr[i], ha='center', va='bottom',fontsize=10)

ax.legend()
fig.suptitle('Farm Power Blockage Comparison',fontsize=20)
# ax.set_title('Farm Power Blockage Comparison',fontsize=20)
fig.tight_layout(rect=(0,0,1,0.93))

plt.show()