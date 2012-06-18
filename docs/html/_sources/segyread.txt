.. documentation for segyread submodule

pygeo.segyread
==============

The :py:mod:`pygeo.segyread` submodule is designed to allow interaction with geophysical (seismic) datafiles that use the `SEG-Y format <http://en.wikipedia.org/wiki/SEG-Y>`_.  The primary purpose of the package is to allow *read-only* access to the SEG-Y data format, though several provisions are made for creating or updating SEG-Y files.

.. automodule:: pygeo.segyread
   :members:

SEGYFile
--------

The  :py:class:`SEGYFile` class represents the SEG-Y or SU datafile efficiently, and initially loads only the metadata necessary to set certain parameters, viz: filesize, endian, data format.  Several objects are created inside the namespace of the :py:class:`SEGYFile` object, viz: **thead**, **bhead**, **trhead**, **endian**, **mendian**, **ns**, **ntr**, **filesize**, **ensembles**.

.. autoclass:: pygeo.segyread.SEGYFile
   :members: __getitem__, findTraces, readTraces, sNormalize, writeFlat, writeSEGY, writeSU

SEGYTraceHeader
---------------

The :py:class:`SEGYTraceHeader` class efficiently indexes the trace headers of the parent :py:class:`SEGYFile`.  This makes it possible to access the headers of an individual trace, or a series of traces without prefetching them from the file on disk.  It interfaces directly with the conventional or memory-mapped file object inside the :py:class:`SEGYFile` object.

.. autoclass:: pygeo.segyread.SEGYTraceHeader
  :members: __getitem__
