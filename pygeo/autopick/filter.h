#include <stdio.h>
#include <omp.h>
#include <Python.h>	
#include <math.h>

void energyRatio (	float inarr[],
			float outarr[],
			Py_ssize_t arrL,
			Py_ssize_t arrW,
			Py_ssize_t strideL,
			Py_ssize_t strideW,
			Py_ssize_t windowsize,
			double damp);
