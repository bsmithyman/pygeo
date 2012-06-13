
// pygeo - a distribution of tools for managing geophysical data
// Copyright (C) 2011, 2012 Brendan Smithyman


// This file is part of pygeo.

// pygeo is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// pygeo is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with pygeo.  If not, see <http://www.gnu.org/licenses/>.

// ----------------------------------------------------------------------

#include "filter.h"

void energyRatio (	float *inarr,
			float *outarr,
			Py_ssize_t arrL,
			Py_ssize_t arrW,
			Py_ssize_t strideL,
			Py_ssize_t strideW,
			Py_ssize_t windowsize,
			double damp) {

  int i, j, pos;
  int jsize = strideL/sizeof(float);
  double totalEnergy, windowTotalEnergy, cursq;

  #pragma omp parallel shared(jsize, inarr, outarr, arrL, arrW, strideL, strideW, windowsize, damp) private(i, j, pos, totalEnergy, windowTotalEnergy, cursq)
  {
    #pragma omp for
    for (i = 0; i < arrL; i++) {
      totalEnergy = 0.;

      // Inner loop over all samples in the trace
      for (j = 0; j < windowsize; j++) {
        pos = i*jsize + j;
        cursq = inarr[pos]*inarr[pos];
        totalEnergy += cursq;
        outarr[pos] = 0.;
      }
      // The total windowTotalEnergy is the totalEnergy right before the window output starts
      windowTotalEnergy = totalEnergy;
      for (j = windowsize; j < arrW; j++) {
        pos = i*jsize + j;
        cursq = inarr[pos]*inarr[pos];
        totalEnergy += cursq;       // Then it begins removing the energy from before the window and adding the energy at the end of the window
        windowTotalEnergy += cursq - inarr[pos - windowsize]*inarr[pos - windowsize];
        // The result is output to the array in parallel
        outarr[pos] = windowTotalEnergy / (totalEnergy + damp);
      } // End local for over j
    } // End parallel for over i
  } // End OMP parallel
}

