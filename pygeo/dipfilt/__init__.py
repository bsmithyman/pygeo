try:
  import pyximport
  pyximport.install()
except:
  print('Cython import failed; pygeo.dipfilt will not function.')
else:
  try:
    from dipfilt import *
  except ImportError:
    print('Could not build/import dipfilt.pyx; pygeo.dipfilt will not function.')
