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
Optimizes wind speed with given wind TI to match FLORIS first row turbines to first row SOWFA powers
"""

input_file="../../example_induction_input.json"
# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file)

# Read in induction options from input file
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction'] = True
fi.IndOpts = Ind_Opts

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with m columns and n rows
m = 3
n = 2

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y])

# Calculate wake
fi.calculate_wake()

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

sowfa_pow_frow = []
for i in range(m):
    sowfa_pow_frow.append(sowfa_pow[i*n])
sowfa_pow_frow = np.array(sowfa_pow_frow)

minval = 5
maxval = 25
x0 = minval
opt_options = {'maxiter': 100,'disp': True,'iprint': 2,'ftol': 1e-7,'eps': 0.01}

def speed_function(x):
    fi.reinitialize_flow_field(wind_speed=x)
    # print(fi.floris.farm.wind_map.input_speed[0])
    fi.calculate_wake()
    turb_powers = np.array(fi.get_turbine_power())/1000

    frow_pow = []
    for i in range(m):
        frow_pow.append(turb_powers[i*n])
    frow_pow = np.array(frow_pow)

    # pow = np.abs(np.sum(np.array(fi.get_turbine_power())/1000 - sowfa_pow))
    return np.linalg.norm(frow_pow-sowfa_pow_frow)/np.linalg.norm(sowfa_pow_frow)*100

plant = minimize(
    speed_function, # the objective function to be minimized
    x0, # Initial guess. Array of real elements the size of n where n in the number of independent variables
    method = 'SLSQP', # SLSQP
    bounds = [(minval,maxval)], # (min,max) pairs for each element in x
    options = opt_options, # maxiter:int and disp (set to true to print convergence messages)
)

print('===================================================')
print('Optimal Wind Speed To Match SOWFA Data')
print('---------------------------------------------------')
print(round(plant.x[0],4))

fi.reinitialize_flow_field(wind_speed=round(plant.x[0],2))
fi.calculate_wake()
turb_powers = np.array(fi.get_turbine_power())/1000

print('===================================================')
print("Farm Power (kW): ", fi.get_farm_power()/1000)
print('===================================================')
print("Turbine Powers:")
print('---------------------------------------------------')
turbine_powers = pd.DataFrame()
for i in range(m):
    temp = []
    for j in range(n):
        temp.append(turb_powers[i*n+j])
    turbine_powers['Column_'+str(i)]=temp
print(turbine_powers)
print('---------------------------------------------------')

print(time.strftime("%H:%M:%S", time.gmtime(time.time()-tstart)))