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


import matplotlib.pyplot as plt

import floris.tools as wfct


# Initialize the FLORIS interface fi
# For basic usage, the florice interface provides a simplified interface to
# the underlying classes
fi = wfct.floris_interface.FlorisInterface("../example_input.json")

# fi.reinitialize_flow_field(layout_array=[[0.0,630.0,1260.0,0.0,630.0,1260.0],[0.0,0.0,0.0,378.0,378.0,378.0]])
# fi.reinitialize_flow_field(layout_array=[[0.0,378.0,756.0,0.0,378.0,756.0],[0.0,0.0,0.0,378.0,378.0,378.0]])
# fi.reinitialize_flow_field(layout_array=[[0.0,882.0,1764.0,0.0,882.0,1764.0],[0.0,0.0,0.0,378.0,378.0,378.0]])

# Calculate wake
fi.calculate_wake()

# Get horizontal plane at default height (hub-height)
hor_plane = fi.get_hor_plane()

# Plot and show
fig, ax = plt.subplots()
wfct.visualization.visualize_cut_plane(hor_plane, ax=ax)
# wfct.visualization.plot_turbines_with_fi(ax,fi)
plt.show()
