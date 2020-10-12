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
import time
# -- Optimization
from scipy.optimize import minimize
from floris.tools.optimization.scipy import optimization

tstart = time.time()

"""
Optimizes turbulence intensity for FLORIS turbine powers to match SOWFA Powers
"""
#Specify SOWFA powers in order of increasing increasing y locations
# sowfa_pow = np.array([1879,793.1,1724.2,767.4,1809.1,798.3]) # 3x2
# sowfa_pow = np.array([1877.1,786.4,784.1,1722,757.3,755.9,1807.2,791.1,787.6]) # 3x3
# sowfa_pow = np.array([1876.9,785,776.8,820,1721.8,755.6,745.6,791.9,1807,789.7,780.5,831.6]) # 3x4
sowfa_pow = np.array([1808.2,778.8,766.7,807.9,817.6,1663.1,720,672,714,779.2,1906.3,802.2,800.6,857.2,900.7]) # 3x5

input_file="../../example_induction_input.json"
# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file)

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with n rows and m columns
n = 5
m = 3

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y], wind_speed= 8.21)

# Read in induction options from input file
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction'] = True
Ind_Opts['Model'] = 'VC'
Ind_Opts['nIter'] = 2
fi.IndOpts = Ind_Opts

# Calculate wake
fi.calculate_wake()

# Store value for total wind farm power
power_initial = fi.get_turbine_power()

print('Initial powers w/ TI = 0.09, U  = 8.19: ', power_initial)

minval = 0.01
maxval = 0.2
x0 = minval
opt_options = {'maxiter': 100,'disp': True,'iprint': 1,'ftol': 1e-8,'eps': 0.01}

def ti_function(x):
    fi.reinitialize_flow_field(turbulence_intensity=x)
    # print(fi.floris.farm.wind_map.input_ti[0])
    fi.calculate_wake()
    # print('Turbine powers with ti of %.2f: ' %x, fi.get_turbine_power())
    pow = np.abs(np.sum(np.array(fi.get_turbine_power())/1000 - sowfa_pow))
    return pow

plant = minimize(
    ti_function, # the objective function to be minimized
    x0, # Initial guess. Array of real elements the size of n where n in the number of independent variables
    method = 'SLSQP', # SLSQP
    bounds = [(minval,maxval) for _ in range(1)], # (min,max) pairs for each element in x
    options = opt_options, # maxiter:int and disp (set to true to print convergence messages)
)

print('===================================================')
print('Optimal TI to Match SOWFA Data')
print('---------------------------------------------------')
print(round(plant.x[0],4))
print('---------------------------------------------------')

# # fi.calculate_wake()
# print(fi.floris.farm.wind_map.input_speed[0])
# print(fi.floris.farm.wind_map.input_ti[0][0])
# print('Turbine powers with %.2f ti: ' %fi.floris.farm.wind_map.input_ti[0], fi.get_turbine_power())

print(time.strftime("%H:%M:%S", time.gmtime(time.time()-tstart)))