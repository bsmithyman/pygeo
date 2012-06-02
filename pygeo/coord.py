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
