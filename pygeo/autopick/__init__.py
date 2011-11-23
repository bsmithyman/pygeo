try:
  import pyximport
  pyximport.install()
except:
  print('Cython import failed; pygeo.autopick will not function.')
else:
  try:
    from autopick import *
  except ImportError:
    print('Coud not build/import autopick.pyx; pygeo.autopick will not function.')
