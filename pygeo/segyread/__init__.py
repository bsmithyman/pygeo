try:
  import pyximport
  pyximport.install()
except:
  print('Cython import failed; pygeo.segyread will use the legacy (pure Python) mode.')
  from segyreadold import *
else:
  try:
    from segyread import *
  except ImportError:
    print('Could not build/import segyread.pyx; pygeo.segyread will use the legacy (pure Python) mode.')
    from segyreadold import *
