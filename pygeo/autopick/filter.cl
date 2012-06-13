
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

typedef struct Params
{
  unsigned int trlen;
  unsigned int jsize;
  unsigned int windowsize;
  float damp;
} Params;

__kernel void energyRatio (	__global float *inarr,
				__global float *outarr,
				__constant struct Params *opts) {

  unsigned int i, j, pos;
  float totalEnergy, windowTotalEnergy, cursq;
  
  i = get_global_id(0);
  totalEnergy = 0.;

  for (j = 0; j < opts->windowsize; j++) { 
    pos = i * opts->jsize + j;
    cursq = inarr[pos]*inarr[pos];
    totalEnergy += cursq;
    outarr[pos] = 0.;
  }

  windowTotalEnergy = totalEnergy;

  for (j = opts->windowsize; j < opts->trlen; j++) {
    pos = i * opts->jsize + j;
    cursq = inarr[pos]*inarr[pos];
    totalEnergy += cursq;
    windowTotalEnergy += cursq - inarr[pos - opts->windowsize]*inarr[pos - opts->windowsize];
    outarr[pos] = windowTotalEnergy / (totalEnergy + opts->damp);
  }
}
