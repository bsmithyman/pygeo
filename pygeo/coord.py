
# pygeo - a distribution of tools for managing geophysical data
# Copyright (C) 2011, 2012 Brendan Smithyman

# This file is part of pygeo.

# pygeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# pygeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with pygeo.  If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------

import numpy as np

def rotationPlan (origpoints, angle=0., about=[0.,0.,0.]):

  # origpoints is an n by 3 array
  newpoints = origpoints.copy()
  temparr = origpoints - np.array(about)
  newpoints[:,0] = np.cos(angle) * temparr[:,0] + np.sin(angle) * temparr[:,1] + about[0]
  newpoints[:,1] = -np.sin(angle) * temparr[:,0] + np.cos(angle) * temparr[:,1] + about[1]
  
  return newpoints

def reduceToLocal (origpoints, angle=0., about=[0.,0.,0.]):
  about = np.array(about)
  return rotationPlan(origpoints, angle, about) - about
