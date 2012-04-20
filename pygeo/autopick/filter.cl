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
