import numpy as np
import pandas as pd
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
import time
# -- Optimization
from scipy.optimize import minimize
from floris.tools.optimization.scipy import optimization

tstart = time.time()

"""
Optimizes heterogeneous input wind speeds to match corresponding first row FLORIS turbine powers to first row 
SOWFA turbine powers using a speedfactor to adjust (slide/translate) upstream 3x1 SOWFA wind velocities and input
into FLORIS and optimizes TI to tune FLORIS to downstream SOWFA powers.
"""

# Read in upstream wind velocities from SOWFA case
InputSpeeds = pd.read_csv('../../../../sowfa_comparisons/SowfaVelocities/SowfaInputVel3x1.csv')
def inSpeeds(numpoints,df):
    locations = np.array([float(i) for i in df[df.NumPoints == numpoints].Locations.values[0].strip('][').split(', ')]) - 320
    velocities = np.array([float(i) for i in df[df.NumPoints == numpoints].Velocities.values[0].strip('][').split(', ')])
    return locations,velocities

# Return wind locations and velocities from sowfa velocity files at sampled points
samploc,sampvel = inSpeeds(19,InputSpeeds)
xwind = np.ones(len(samploc))*(-252)

input_file="../../../example_induction_input.json"
# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file)

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with m columns and n rows
m = 3
n = 5

#Specify SOWFA powers in order of increasing increasing y locations
if m == 3 and n == 1:
    sowfa_pow = np.array([1888.4,1736.6,1818.6]) # 3x1
if m == 3 and n == 2:
    sowfa_pow = np.array([1879,793.1,1724.2,767.4,1809.1,798.3]) # 3x2
if m == 3 and n == 3:
    sowfa_pow = np.array([1877.1,786.4,784.1,1722,757.3,755.9,1807.2,791.1,787.6]) # 3x3
if m == 3 and n == 4:
    sowfa_pow = np.array([1876.9,785,776.8,820,1721.8,755.6,745.6,791.9,1807,789.7,780.5,831.6]) # 3x4
if m == 3 and n == 5:
    sowfa_pow = np.array([1808.2,778.8,766.7,807.9,817.6,1663.1,720,672,714,779.2,1906.3,802.2,800.6,857.2,900.7]) # 3x5

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y], wind_speed=list(sampvel),wind_layout=[xwind,list(samploc)])

# Read in induction options from input file
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction'] = True
fi.IndOpts = Ind_Opts

ti_multfact = 100 # Multiply ti parameter to be on same order as windspeed
ti_minval = 0.01 * ti_multfact
ti_maxval = 0.2 * ti_multfact
minspeedfac = -3
maxspeedfac = 3
x0 = [ti_minval, minspeedfac]
opt_options = {'maxiter': 200,'disp': True,'iprint': 2,'ftol': 1e-8,'eps': 0.01}

def opt_function(param):
    fi.reinitialize_flow_field(turbulence_intensity=param[0]/ti_multfact,wind_speed=list(sampvel+param[1]))
    fi.calculate_wake()
    turb_powers = np.array(fi.get_turbine_power())/1000
    # pow = np.sum(np.abs((turb_powers-sowfa_pow)/sowfa_pow))
    pow = np.linalg.norm(turb_powers - sowfa_pow)/np.linalg.norm(sowfa_pow)
    return pow

plant = minimize(
    opt_function, # the objective function to be minimized
    x0, # Initial guess. Array of real elements the size of n where n in the number of independent variables
    method = 'SLSQP', # SLSQP
    bounds = [(ti_minval,ti_maxval),(minspeedfac,maxspeedfac)], # (min,max) pairs for each element in x
    options = opt_options, # maxiter:int and disp (set to true to print convergence messages)
)

# print('===================================================')
# print('Optimal TI to Match SOWFA Data')
# print('---------------------------------------------------')
# # print(plant.x[0]/ti_multfact)
# print('Rounded TI: ', round(plant.x[0]/ti_multfact,2))

# Reinitialize flow with optimized turbulence intensity and calculate turbine powers
fi.reinitialize_flow_field(turbulence_intensity=round(plant.x[0]/ti_multfact,4),wind_speed=list(sampvel+plant.x[1]))
fi.calculate_wake()

# hor_plane = fi.get_hor_plane()
# # Plot 
# fig,ax = plt.subplots()
# wfct.visualization.visualize_cut_plane(hor_plane,ax=ax,fig=fig,cbar=True)
# wfct.visualization.plot_turbines_with_fi(ax,fi)

print('===================================================')
print('Optimal Speed Factor: ', round(plant.x[1],4))
print('==============================================')
print('Turbine ti: ', fi.floris.farm.wind_map.input_ti[0])
print('==============================================')
print("Farm Power (kW): ", fi.get_farm_power()/1000)
print('==============================================')
print("Turbine Powers:")
print('----------------------------------------------')
turbpow = fi.get_turbine_power()
turbine_powers = pd.DataFrame()
for i in range(m):
    temp = []
    for j in range(n):
        temp.append(turbpow[i*n+j]/1000)
    turbine_powers['Column_'+str(i)]=temp
print(turbine_powers)

# print(turbpow)

print(time.strftime("%H:%M:%S", time.gmtime(time.time()-tstart)))
plt.show()