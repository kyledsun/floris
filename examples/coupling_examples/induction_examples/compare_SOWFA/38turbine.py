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
Comparison of Baseline and Blockage FLORIS to 38 turbine SOWFA Case.
"""

horizontalPlane = False

input_file="../../example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# case1_baseline
layout_x = list((2436.8, 3082.4, 2759.6, 2114.0, 1791.2, 2114.0, 2759.6, 3728.0, 3555.0, 3082.4, 2436.8, 1791.2, 1318.6, 1145.6, 1318.6, 1791.2, 2436.8, 3082.4, 3555.0, 4373.6, 4268.7, 3965.2, 3496.1, 2912.3, 2276.9, 1658.8, 1125.0, 733.4, 526.4, 526.4, 733.4, 1125.0, 1658.8, 2276.9, 2912.3, 3496.1, 3965.2, 4268.7))
layout_y = list((2436.8, 2436.8, 2995.9, 2995.0, 2436.8, 1877.0, 1877.7, 2436.0, 3082.4, 3555.0, 3728.0, 3555.0, 3082.4, 2436.8, 1791.2, 1318.6, 1145.6, 1318.6, 1791.2, 2436.8, 3065.7, 3626.4, 4058.2, 4314.3, 4367.0, 4210.5, 3861.8, 3358.6, 2755.6, 2118.0, 1515.0, 1011.8, 663.1, 506.6, 559.3, 815.4, 1247.2, 1807.9))
sowfa_pow = np.array([782.2,879.8,606.1,1550,626.7,1594.6,622.3,826.8,1033.8,1052.1,1936.4,1971.7,2023.2,1952.7,1813.4,1954.8,1818.4,1078.3,845.9,896.9,685.6,924.1,1907.6,822.2,1971,1871.3,2016.4,2089.6,2046.3,1966.7,1779.9,1939.4,1988.3,1756.3,782.8,1789.6,1000.9,636.9])

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y],wind_speed=8.39,turbulence_intensity=0.065)

# Read in induction options
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction']=True
Ind_Opts["Model"] = 'VC'
fi.IndOpts = Ind_Opts

# Calculate wake
fi.calculate_wake()
#------------------------------------------- Baseline Power -------------------------------------------
fi_base = copy.deepcopy(fi)

# Read in induction options
Ind_Opts = fi_base.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction']=False
fi_base.IndOpts = Ind_Opts

# Calculate wake
fi_base.calculate_wake()

print('==============================================')
sowfaFarmPow = np.sum(sowfa_pow)
print("SOWFA Farm Power (kW): %.2f" %sowfaFarmPow)
print("Baseline Farm Power (kW): %.2f" %(fi_base.get_farm_power()/1000))
print("Baseline Farm Power Error (%%): %.2f%%" %((fi_base.get_farm_power()/1000-sowfaFarmPow)/sowfaFarmPow*100))
print('----------------------------------------------')
print("Blockage Farm Power (kW): %.2f" %(fi.get_farm_power()/1000))
print("Blockage Farm Power Error (%%): %.2f%%" %((fi.get_farm_power()/1000-sowfaFarmPow)/sowfaFarmPow*100))
print('==============================================')

# Turbine Powers
baseturbpow = np.array(fi_base.get_turbine_power())/1000 # Returns turbine powers and converts to kW
blockturbpow = np.array(fi.get_turbine_power())/1000 # Returns turbine powers and converts to kW

# Get tubine locations
df = pd.DataFrame(columns=['xloc','yloc'])
xloc,yloc = fi.get_turbine_layout()
df['xloc'] = xloc
df['yloc'] = yloc
df['SOWFA'] = sowfa_pow
df['FLORIS_baseline'] = baseturbpow
df['FLORIS_baseline_Err'] = (baseturbpow-sowfa_pow)/sowfa_pow*100
df['FLORIS_blockage'] = blockturbpow
df['FLORIS_blockage_Err'] = (blockturbpow-sowfa_pow)/sowfa_pow*100
df['BW'] = ['worse' for i in range(df.shape[0])]
df.at[df[abs(df.FLORIS_baseline_Err) > abs(df.FLORIS_blockage_Err)].index,'BW'] = 'better'

# print(df)

fig,ax=plt.subplots(1,2,figsize=(10,5))
im = ax[0].scatter(df.xloc,df.yloc,c=df.FLORIS_baseline_Err,cmap='RdBu',vmin=-30,vmax=30)
ax[0].set_title('FLORIS Baseline--SOWFA Error',fontsize=14)
divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right",size='2%',pad = 0.07)
cb = fig.colorbar(im, cax=cax)
cax.set_ylabel('Percent Error (%)',rotation=270,fontsize=12,labelpad=18)
ax[0].set_aspect('equal')

im = ax[1].scatter(df.xloc,df.yloc,c=df.FLORIS_blockage_Err,cmap='RdBu',vmin=-30,vmax=30)
# im = ax[1].scatter(df.xloc,df.yloc,c=df.FLORIS_blockage_Err,cmap='nipy_spectral',vmin=-30,vmax=30)
ax[1].set_title('FLORIS Blockage--SOWFA Error',fontsize=14)
divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right",size='2%',pad = 0.07)
cb = fig.colorbar(im, cax=cax)
cax.set_ylabel('Percent Error (%)',rotation=270,fontsize=12,labelpad=18)
ax[1].set_aspect('equal')
fig.tight_layout()

# fig,ax = plt.subplots()
# ax.scatter(df[df.BW=='worse'].xloc,df[df.BW=='worse'].yloc,label='Baseline')
# ax.scatter(df[df.BW=='better'].xloc,df[df.BW=='better'].yloc,label='Blockage')
# ax.legend()
# ax.set_aspect('equal')
# ax.set_title('Better Prediction Compared to SOWFA',fontsize=16)


if horizontalPlane:
    # Initialize the horizontal cut
    hor_plane = fi_base.get_hor_plane(x_bounds=[-1000,6000],x_resolution=300)
    # Initialize the horizontal cut
    hor_plane_block = fi.get_hor_plane(x_bounds=[-1000,6000],x_resolution=300)

    # Plot and Show
    fig, ax = plt.subplots(1,2)
    # wfct.visualization.visualize_cut_plane(hor_plane, ax=ax, fig=fig, cbar=True,levels=levelsLines,cmap='gist_stern_r',minSpeed=0,maxSpeed=8.0)
    wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0], fig=fig, cbar=True,levels=None)#,minSpeed=0,maxSpeed=8.0)
    wfct.visualization.visualize_cut_plane(hor_plane_block, ax=ax[1], fig=fig, cbar=True,levels=np.arange(0.97,1.01,0.005) * fi.floris.farm.wind_map.input_speed[0])#,minSpeed=0,maxSpeed=8.0)
    wfct.visualization.plot_turbines_with_fi(ax[0],fi_base)
    wfct.visualization.plot_turbines_with_fi(ax[1],fi)
    ax[0].set_title('Baseline Floris',fontsize=14)
    ax[1].set_title('Floris with Blockage',fontsize=14)

print('----------------------------------------------')
print('Time Elapsed: ', (time.time()-tstart))
plt.show()