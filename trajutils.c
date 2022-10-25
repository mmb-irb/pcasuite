#include "pcazip.h"
#include "trajutils.h"

/* This method receives a trajectory (x), a reference frame (ref) and
 * modifies the trajectory to the RMSd best fit with the reference frame
 * using the Kabsch algorithm (standard version) */
void trajfit (float *x, float *ref, int natoms, int nframes) {
  int findmove, i, j;
  float xt, yt, zt, rmsd;
  float r[9], v[3];

  findmove = TRUE;
  
  for (i = 0; i < nframes; i++) {
    matfit(&natoms, ref, &x[i*natoms*3], r, v, &rmsd, &findmove);
    for (j = 0; j < natoms; j++) {
      xt = r[0] * x[0+((j+(i*natoms))*3)] +
           r[3] * x[1+((j+(i*natoms))*3)] +
           r[6] * x[2+((j+(i*natoms))*3)] + v[0];
      yt = r[1] * x[0+((j+(i*natoms))*3)] +
           r[4] * x[1+((j+(i*natoms))*3)] +
           r[7] * x[2+((j+(i*natoms))*3)] + v[1];
      zt = r[2] * x[0+((j+(i*natoms))*3)] +
           r[5] * x[1+((j+(i*natoms))*3)] +
           r[8] * x[2+((j+(i*natoms))*3)] + v[2];
      x[0+((j+(i*natoms))*3)] = xt;
      x[1+((j+(i*natoms))*3)] = yt;
      x[2+((j+(i*natoms))*3)] = zt;
    }
  }
}

/* This method receives a trajectory (x) and computes the mean structure,
 * leaving the computed mean structure in avg parameter */
void trajavg (float *x, float *avg, int natoms, int nframes) {
  int i, j;

  /* Average vector initialization */
  for (j = 0; j < natoms * 3; j++) {
    avg[j] = 0.0;
  }

  /* Sum of data in different frames */
  for (i = 0; i < nframes; i++) {
    for (j = 0; j < natoms * 3; j++ ) {
      avg[j] += x[((natoms * 3) * i) + j];
    }
  }

  /* Mean division */
  for (j = 0; j < natoms * 3; j++) {
    avg[j] /= nframes;
  }
}
