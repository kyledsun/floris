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
Optimizes heterogeneous input wind speed (at 3 points specified upstream of each turbine column)
to match first row turbine powers between FLORIS and SOWFA, then optimizes TI to match FLORIS 
turbine powers throughout the farm to SOWFA turbine powers.
"""

input_file="../../../example_induction_input.json"
# Initialize floris interface object
fi = wfct.floris_interface.FlorisInterface(input_file)

# Read in induction options from input file
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
Ind_Opts['induction'] = True
Ind_Opts['Model'] = 'VC'
Ind_Opts['nIter'] = 2
fi.IndOpts = Ind_Opts

# Read in input velocities from 3x1 SOWFA case
InputSpeeds = pd.read_csv('../../../../sowfa_comparisons/SowfaVelocities/SowfaInputVel3x1.csv')
def inSpeeds(numpoints,df):
    locations = np.array([float(i) for i in df[df.NumPoints == numpoints].Locations.values[0].strip('][').split(', ')]) - 320
    velocities = np.array([float(i) for i in df[df.NumPoints == numpoints].Velocities.values[0].strip('][').split(', ')]) #* (8.38 / 8)
    return locations,velocities

# Return wind locations and velocities from sowfa velocity files at sampled points
samploc,sampvel = inSpeeds(3,InputSpeeds)
xwind = np.ones(len(samploc))*-252

# Set paramters for iteration test
sep = 5 # streamwise separation for turbines (*D)
sepy = 3 # spanwise spearation between turbines (*D)
# Creates a turbine field with m columns and n rows
m = 3
n = 1

D = fi.floris.farm.turbines[0].rotor_diameter
layout_x = []
layout_y = []
for i in range(m):
    for j in range(n):
        layout_x.append(j*sep*D)
        layout_y.append(i*sepy*D)

# Reinitialize flow field with new specified layout
fi.reinitialize_flow_field(layout_array=[layout_x,layout_y], wind_speed=list(sampvel),wind_layout=[xwind,list(samploc)])
# Calculate wake
fi.calculate_wake()

#Specify SOWFA powers in order of increasing increasing y,x locations
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

# Collect Turbine Powers for First Row Turbines
sowfa_pow_frow = []
for i in range(m):
    sowfa_pow_frow.append(sowfa_pow[i*n])
sowfa_pow_frow = np.array(sowfa_pow_frow)

minspeedfac = -3
maxspeedfac = 3
x0 = [minspeedfac,minspeedfac,minspeedfac]
opt_options = {'maxiter': 100,'disp': True,'iprint': 2,'ftol': 1e-7,'eps': 0.01}
# opt_options = {'maxiter': 100,'disp': True,'iprint': 2,'ftol': 1e-7,'eps': 0.05}

# Optimization function to optimize each of three heterogeneous input wind speeds individually.
def speed_func(param):
    fi.reinitialize_flow_field(wind_speed=[sampvel[0]+param[0],sampvel[1]+param[1],sampvel[2]+param[2]])
    fi.calculate_wake()
    turb_powers = np.array(fi.get_turbine_power())/1000

    frow_pow = []
    for i in range(m):
        frow_pow.append(turb_powers[i*n])
    frow_pow = np.array(frow_pow)

    return np.linalg.norm(frow_pow-sowfa_pow_frow)/np.linalg.norm(sowfa_pow_frow)*100

plant = minimize(
    speed_func, # the objective function to be minimized
    x0, # Initial guess. Array of real elements the size of n where n in the number of independent variables
    method = 'SLSQP', # SLSQP
    bounds = [(minspeedfac,maxspeedfac),(minspeedfac,maxspeedfac),(minspeedfac,maxspeedfac)], # (min,max) pairs for each element in x
    options = opt_options, # maxiter:int and disp (set to true to print convergence messages)
)

fi.reinitialize_flow_field(wind_speed=[sampvel[0]+round(plant.x[0],4),sampvel[1]+round(plant.x[1],4),sampvel[2]+round(plant.x[2],4)])
fi.calculate_wake()
turb_powers = np.array(fi.get_turbine_power())/1000

multfact = 100
minval = 0.06 * multfact
maxval = 0.2 * multfact
x0 = minval
opt_options = {'maxiter': 100,'disp': True,'iprint': 2,'ftol': 1e-8,'eps': 0.01}

# Optimization function to optimize TI to tune FLORIS to match SOWFA turbine powers
def ti_function(x):
    fi.reinitialize_flow_field(turbulence_intensity=x/multfact)
    fi.calculate_wake()
    turb_powers = np.array(fi.get_turbine_power())/1000
    # pow = np.abs(np.sum(np.array(fi.get_turbine_power())/1000 - sowfa_pow))
    return np.linalg.norm(turb_powers-sowfa_pow)/np.linalg.norm(sowfa_pow)

plant2 = minimize(
    ti_function, # the objective function to be minimized
    x0, # Initial guess. Array of real elements the size of n where n in the number of independent variables
    method = 'SLSQP', # SLSQP
    bounds = [(minval,maxval)], # (min,max) pairs for each element in x
    options = opt_options, # maxiter:int and disp (set to true to print convergence messages)
)

fi.reinitialize_flow_field(wind_speed=[sampvel[0]+round(plant.x[0],4),sampvel[1]+round(plant.x[1],4),sampvel[2]+round(plant.x[2],4)],turbulence_intensity=round(plant2.x[0]/multfact,3))
fi.calculate_wake()
turb_powers = np.array(fi.get_turbine_power())/1000

# hor_plane = fi.get_hor_plane()
# fig,ax = plt.subplots()
# wfct.visualization.visualize_cut_plane(hor_plane,ax=ax,fig=fig,cbar=True)
# wfct.visualization.plot_turbines_with_fi(ax,fi)
# plt.show()

print('===================================================')
print('Optimal Wind Speed Factors To Match SOWFA Data')
print('---------------------------------------------------')
speedfactors = pd.DataFrame([[round(plant.x[0],4),round(plant.x[1],4),round(plant.x[2],4)]],index=['speedfactor','windspeeds'],columns=['Column 1','Column 2','Column 3'])
speedfactors.iloc[1] = [sampvel[0]+round(plant.x[0],4),sampvel[1]+round(plant.x[1],4),sampvel[2]+round(plant.x[2],4)]
print(speedfactors)
print('===================================================')
print('Optimal TI Factor To Match SOWFA Data')
print('---------------------------------------------------')
print(round(plant2.x[0]/multfact,3))
print('===================================================')
print("Farm Power: ", fi.get_farm_power())
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