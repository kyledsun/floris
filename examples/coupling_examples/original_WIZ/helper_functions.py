import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable


# --- Local
import floris.tools as wfct
from floris.utilities import Vec3

def get_HH_plane_vel(input_file, input_dict, Induction=False, resolution=Vec3(232, 114, 3),bounds_to_set=None):
    fi = wfct.floris_interface.FlorisInterface(input_file, input_dict)
    # print("Ind Opts:", Ind_Opts)
    Ind_Opts = fi.floris.farm.flow_field.Ind_Opts
    if Induction:
        print(">>> Blockage effect is modeled")
        Ind_Opts['induction'] = True
    else:
        Ind_Opts['induction'] = False

    # fi.reinitialize_flow_field(layout_array=[[0],[0]])
    # Calculate wake
    # fi.calculate_wake(yaw_angles=[15],Ind_Opts=Ind_Opts)
    #fi.calculate_wake(yaw_angles=[15,5,0.0,15,5,0],Ind_Opts=Ind_Opts)
    #fd=fi.get_flow_data(resolution = resolution, Ind_Opts=Ind_Opts)

    # # Initialize the horizontal cut
    hor_plane = fi.get_hor_plane(
        height = fi.floris.farm.turbines[0].hub_height,
        x_resolution = resolution.x1,
        y_resolution = resolution.x2,
        x_bounds = tuple(bounds_to_set[0:2]),
        y_bounds = tuple(bounds_to_set[2:4]),
        Ind_Opts = Ind_Opts
        )
    u_mesh = hor_plane.df.u.values.reshape(resolution.x2,resolution.x1)
    v_mesh = hor_plane.df.v.values.reshape(resolution.x2,resolution.x1)
    
    return np.unique(hor_plane.df.x1),np.unique(hor_plane.df.x2), u_mesh, v_mesh, fi


def savePlane(x,y=None,u=None,v=None,base='',U0=None):
    if u is None:
        p=x
        base=y
        x=p[0]
        y=p[1]
        u=p[2]
        v=p[3]
    if U0 is not None:
        u=u/U0
        v=v/U0

    base = os.path.normpath(base)
    np.save(base+'x.npy',x)
    np.save(base+'y.npy',y)
    np.save(base+'u.npy',u)
    np.save(base+'v.npy',v)


def loadPlane(base=''):
    base = os.path.normpath(base)
    x = np.load(base+'x.npy')
    y = np.load(base+'y.npy')
    u = np.load(base+'u.npy')
    v = np.load(base+'v.npy')
    return x,y,u,v


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
#         qv=streamQuiver(ax,sp,n=5,scale=40,angles='xy')

    ax.set_aspect('equal')

    return im
