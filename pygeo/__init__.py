try:
  import pyximport
  pyximport.install()
except:
  print('Cython import failed; the dipfilt sub-module will not function.')
else:
  import dipfilt
