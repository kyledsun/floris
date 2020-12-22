# Copyright 2020 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# See https://floris.readthedocs.io for documentation

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import floris.tools as wfct

"""
Visualizes horizontal plane for specified SOWFA simulation, plots most 
upstream wind velocities at turbine hub height and samples the upstream
velocities to be used in FLORIS heterogeneous input comparisons.
"""

sampledPoints = False
savesamples = False
sowfa_case = wfct.sowfa_utilities.SowfaInterface("../../../../../for_chris/for_chris/case_7_3_turbine") # 3x1
# sowfa_case = wfct.sowfa_utilities.SowfaInterface("../../../../../for_chris/for_chris/case_5_6_turbine") # 3x2
# sowfa_case = wfct.sowfa_utilities.SowfaInterface("../../../../../for_chris/for_chris/case_3_9_turbine") # 3x3
# sowfa_case = wfct.sowfa_utilities.SowfaInterface("../../../../../for_chris/for_chris/case_1_12_turbine") # 3x4
# sowfa_case = wfct.sowfa_utilities.SowfaInterface("../../../../../for_chris/for_chris/case_1_baseline") # 3x5

# Demonstrate flow field visualizations

# Get the horizontal cut plane at 90 m hub height
hor_plane = sowfa_case.get_hor_plane(90)
# Get the cross plane cut at x = 0 (input conditions)
cross_plane = sowfa_case.get_cross_plane(0)

# Plot cut planes
fig,ax = plt.subplots(figsize=(7,2.75))
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax, fig=fig, cbar=True,minSpeed=4,maxSpeed=9)
ax.set_title("SOWFA Streamwise Velocity at Hub Height", fontsize=20)

# fig, ax = plt.subplots(nrows = 2)
# wfct.visualization.visualize_cut_plane(hor_plane, ax=ax[0], fig=fig, cbar=True)
# ax[0].set_title("Horizontal Plane at Hub Height")
# wfct.visualization.visualize_cut_plane(cross_plane, ax=ax[1])
# ax[1].set_title("Cross Plane at x = 0")

# Save y locations of turbine columns
y_locations = np.unique(sowfa_case.layout_y)
x_locations = np.unique(sowfa_case.layout_x)
print('x_locations: ', x_locations)
# Save line of crossplane at hub height
cplane = cross_plane.df[cross_plane.df.x2 == 90.0]

# Plot streamwise velocities at hub height at input
fig,ax = plt.subplots()
ax.scatter(cplane.x1, cplane.u, s=5)
ax.plot(np.ones(2)*y_locations[0],[cplane.u.min(),cplane.u.max()])
ax.plot(np.ones(2)*y_locations[1],[cplane.u.min(),cplane.u.max()])
ax.plot(np.ones(2)*y_locations[2],[cplane.u.min(),cplane.u.max()])
ax.set_xlabel('y [-]')
ax.set_ylabel('stream wise velocity (u)')
ax.set_title('Wind Speed at Hub Height at x=0')

# Function to find nearest value in an array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

print('-----------------------------------------------')
print('y_locations: ', y_locations)
print('-----------------------------------------------')
print('Velocities upstream of Turbines')
print('-----------------------------------------------')
for i in range(len(y_locations)):
    yloc = find_nearest(cplane.x1,y_locations[i])
    print('Turbine Column at y at %.1f = ' %yloc, round(cplane[cplane.x1==yloc].u.values[0],4))

def findpoints(locations,df):
    locs=[]; samplepts=[]
    for i in locations:
        temploc = find_nearest(df.x1,i)
        locs.append(temploc)
        samplepts.append(df[df.x1==temploc].u.values[0])
    return locs,samplepts

print('Max: ',cplane.x1.max())
print('Min: ',cplane.x1.min())

# Print sampled points to specify flow in FLORIS
if sampledPoints:
    fig,ax = plt.subplots()
    # All Points
    ax.plot(cplane.x1, cplane.u, c='black',label='SOWFA Velocity',zorder=1)
    # 3 Sampled Points
    yloc3,sampvel3=findpoints(y_locations,cplane)
    ax.scatter(yloc3, sampvel3, c='C2', label = '3 points',zorder=2, s=125)
    # 9 Sampled Points
    yloc9,sampvel9=findpoints([0,170,320,520,700,870,1080,1220,1400],cplane)
    # yloc9,sampvel9=findpoints([0,210,520,660,900,1110,1280,1560,1780],cplane) #3x5 case
    ax.scatter(yloc9, sampvel9, c='C1', label = '9 points',zorder=3,s=75)
    # 19 Sampled Points
    yloc19,sampvel19=findpoints([0,80,160,230,320,390,470,540,620,700,780,860,930,1010,1080,1170,1240,1320,1400],cplane)
    # yloc19,sampvel19=findpoints([0,100,200,300,400,520,600,690,790,900,990,1100,1190,1280,1380,1480,1580,1680,1780],cplane) # 3x5 case
    ax.scatter(yloc19, sampvel19, c='C3', label = '19 points',zorder=4,s=25)

    ax.legend()
    ax.set_xlabel('y [-]')
    ax.set_ylabel('stream wise velocity (u)')
    ax.set_title('Wind Speed at Hub Height at x=0')

    # Save sampled points into csv file
    if savesamples:
        df = pd.DataFrame(columns = ['NumPoints','Locations','Velocities'])
        df['NumPoints'] = [3,9,19,100]
        df['Locations'] = [yloc3,yloc9,yloc19,[i for i in cplane.x1]]
        df['Velocities'] = [sampvel3,sampvel9,sampvel19,[i for i in cplane.u]]
        df.to_csv('SowfaVelocities/SowfaInputVel3x5.csv',index=None)

plt.show()