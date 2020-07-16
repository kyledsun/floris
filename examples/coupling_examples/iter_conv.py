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

fontsize = 12
plt.rc('font', family='serif')
plt.rc('font', size=12)

# --- Plot options
nStreamlines=0
U0 = 8
minSpeed=0.5
maxSpeed=1.03
levelsLines=np.sort([1.05,1.0,0.99,0.98,0.95,0.9,0.5])
horPlot = False

# --- Resolution Parameters
ny=30 # 200
nx=ny*4

input_file="OptLayout_2x3.json"
input_dict=None
# D=126
resolution=Vec3(nx, ny, 2)
# bounds=[-4*D,16*D,-2*D-10,5*D+10,89,90] # xmin xmax .. zmin zmax
Induction = True

# Helper function to plot wind farm velocity fields
def plotPlane(x,y,u,v,ax,minSpeed=None,maxSpeed=None, cmap='coolwarm', colors='w', linewidths=0.8, alpha=0.3,nStreamlines=0,levelsContour=None,levelsLines=None,axial=False):
    if axial:
        Speed=u
    else:
        Speed=np.sqrt(u**2+v**2)
    if minSpeed is None:
        minSpeed = Speed.min()
        maxSpeed = Speed.max()

    if not ax:
        fig, ax = plt.subplots()
    Z = np.ma.masked_where(np.isnan(Speed), Speed)

    # Plot the cut-through
    im = ax.pcolormesh(x, y, Z, cmap=cmap, vmin=minSpeed, vmax=maxSpeed)
    if levelsLines is not None:
        rcParams['contour.negative_linestyle'] = 'solid'
        cs=ax.contour(x, y, Z, levels=levelsLines, colors=colors, linewidths=linewidths, alpha=alpha)

    if nStreamlines>0:
        yseed=np.linspace(min(y)*0.9,max(y)*0.9,nStreamlines)
        start=np.array([yseed*0,yseed])
        sp=ax.streamplot(x,y,u,v,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
    ax.set_aspect('equal')

    return im

# Initialize the floris interface
fi = wfct.floris_interface.FlorisInterface(input_file, input_dict)

# Read in induction options from input file
Ind_Opts = fi.floris.farm.flow_field.Ind_Opts

if Induction:
    print(">>> Blockage effect is modeled")
    Ind_Opts['induction'] = True
else:
    Ind_Opts['induction'] = False

# Turbine parameters
D = fi.floris.farm.turbines[0].rotor_diameter
U0 = fi.floris.farm.initwind[0]

layout_y = [0,0,0,0,0]

# Array of diameter separations between turbines
sep = [3,5,7]
titles = ['3D Seperation','5D Separation','7D Separation']

# Number of iterations used
nIter = np.arange(1,6)

print('Number of iterations:', Ind_Opts['nIter'])

sep3D = pd.DataFrame(index=nIter,columns = ['turbine '+str(i) for i in range(len(layout_y))])
sep5D = pd.DataFrame(index=nIter,columns = ['turbine '+str(i) for i in range(len(layout_y))])
sep7D = pd.DataFrame(index=nIter,columns = ['turbine '+str(i) for i in range(len(layout_y))])

for j in range(len(sep)):
    layout_x = [0,sep[j]*D,2*sep[j]*D,3*sep[j]*D,4*sep[j]*D]
    bounds=[-4*D,sep[j]*5.5*D,-2*D-10,2*D+10,89,90] # xmin xmax .. zmin zmax
    for i in range(len(nIter)):
        Ind_Opts['nIter'] = nIter[i]
        # print('Number of iterations:', Ind_Opts['nIter'])

        # Reinitialize flow flow field with specified layout
        fi.reinitialize_flow_field(layout_array=[layout_x, layout_y])

        # Calculate wake and get horizontal plane at turbine height for original farm field
        fi.calculate_wake(Ind_Opts=Ind_Opts)

        # Store value for total wind farm power
        power_initial = fi.get_farm_power()
        turbine_powers = fi.get_turbine_power()

        if j == 0:
            sep3D.iloc[i] = turbine_powers
        elif j == 1:
            sep5D.iloc[i] = turbine_powers
        elif j == 2:
            sep7D.iloc[i] = turbine_powers
    
    # Plot horizontal cut of the wind farm for final iteration for each turbine layout
    if horPlot:
        hor_plane = fi.get_hor_plane(
                height = fi.floris.farm.turbines[0].hub_height,
                x_resolution = resolution.x1,
                y_resolution = resolution.x2,
                x_bounds = tuple(bounds[0:2]),
                y_bounds = tuple(bounds[2:4]),
                Ind_Opts = Ind_Opts
                )
        u_mesh = hor_plane.df.u.values.reshape(resolution.x2,resolution.x1)
        v_mesh = hor_plane.df.v.values.reshape(resolution.x2,resolution.x1)

        if U0 is not None:
                u_mesh=u_mesh/U0
                v_mesh=v_mesh/U0

        plane_x, plane_y = np.unique(hor_plane.df.x1),np.unique(hor_plane.df.x2)

        # --- Plot and show
        fig, ax = plt.subplots()
        im = plotPlane(plane_x/D,plane_y/D,u_mesh,v_mesh,ax,minSpeed=minSpeed,maxSpeed=maxSpeed,
                    nStreamlines=nStreamlines,levelsLines=levelsLines,axial=True,colors ='k')
        ax.title.set_text(titles[j])
        wfct.visualization.plot_turbines_with_fi(ax,fi)
        ax.set_ylabel('r/D [-]')
        ax.set_xlim([-4,sep[j]*5.5])
        ax.set_ylim([-2,2])

# # Plots the powers of turbines 1 through 4 for 3D separation as a function of the number of iterations
# fig,ax = plt.subplots()
# ax.plot(sep3D.reindex(columns = ['turbine '+str(i) for i in range(1,len(layout_y))]))
# ax.set_ylabel('Power')
# ax.set_xlabel('Number of Iterations')
# ax.set_title('3D Separation')
# ax.legend(['turbine '+str(i) for i in range(1,len(layout_y))])

# Plots the power for each turbine individually for 5D separation as a function of nIter
for i in range(len(layout_y)):
    fig,ax = plt.subplots()
    ax.plot(sep3D['turbine '+str(i)])
    ax.set_title('Turbine '+str(i))
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Power')

# # Plots the power of each turbine as a function of separation and iterations
# for i in range(len(layout_y)):
#     fig,ax = plt.subplots()
#     ax.plot(sep3D['turbine '+str(i)])
#     ax.plot(sep5D['turbine '+str(i)])
#     ax.plot(sep7D['turbine '+str(i)])
#     ax.set_title('Turbine '+str(i))
#     ax.set_xlabel('Iterations')
#     ax.set_ylabel('Power')
#     ax.legend(['3D','5D','7D'])
#     plt.yscale('log')

print('3D Seperation Between Turbines:\n',sep3D)
print('5D Seperation Between Turbines:\n',sep5D)
print('7D Seperation Between Turbines:\n',sep7D)

plt.show()