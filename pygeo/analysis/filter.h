
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

#include <omp.h>
#include <Python.h>	
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void energyRatio (	float inarr[],
			float outarr[],
			Py_ssize_t arrL,
			Py_ssize_t arrW,
			Py_ssize_t strideL,
			Py_ssize_t strideW,
			Py_ssize_t windowsize,
			double damp);

void traceNormalize (	float inarr[],
			float outarr[],
			Py_ssize_t arrL,
			Py_ssize_t arrW,
			Py_ssize_t strideL,
			Py_ssize_t strideW);
