import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from floris.induction import options_dict
from helper_functions import *

fontsize = 12
plt.rc('font', family='serif')
plt.rc('font', size=12)

# --- Plot options
bExport=False
bCompute=True
bColorBar=False
nStreamlines=0
U0 =6
minSpeed=0.5
maxSpeed=1.03
levelsLines=np.sort([1.05,1.0,0.99,0.98,0.95,0.9,0.5])

# ---
ny=30 # 200
nx=ny*4

input_file="Layout_1x2.json"
input_dict=None
D=126
resolution=Vec3(nx, ny, 2)
bounds=[-4*D,14*D,-2*D-10,2*D+10,89,90] # xmin xmax .. zmin zmax
Ind_Opts=options_dict()
Ind_Opts['Rfact']=1.0
Ind_Opts['GammaFact']=1.0
Ind_Opts['Ground']=True

fi = wfct.floris_interface.FlorisInterface(input_file, input_dict)
# Calculate wake
fi.calculate_wake(Ind_Opts=Ind_Opts)

bounds_to_set = bounds

hor_plane = fi.get_hor_plane(
    height = fi.floris.farm.turbines[0].hub_height,
    x_resolution = resolution.x1,
    y_resolution = resolution.x2,
    x_bounds = tuple(bounds_to_set[0:2]),
    y_bounds = tuple(bounds_to_set[2:4])
    )
u_mesh = hor_plane.df.u.values.reshape(resolution.x2,resolution.x1)
v_mesh = hor_plane.df.v.values.reshape(resolution.x2,resolution.x1)

planes.apppend(np.unique(hor_plane.df.x1),np.unique(hor_plane.df.x2), u_mesh, v_mesh, fi)