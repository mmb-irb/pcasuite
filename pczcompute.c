#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "x_io.h"
#include "pczcompute.h"
#include "trajutils.h"
#include "utils.h"
#include "constants.h"

/* This method computes the atomic fluctuations for a given eigenvector,
 * or for all the eigenvectors, if none (-1) is given */
float *pczfluc (PCZFILE *pczf, int iv) {
  float *fluc, *evec, *avg, *snapshot;
  float x, y, z;
  int i, j;

  fluc = (float *)smalloc (pczf->natoms * sizeof (float));
  if (iv >= 0) {
    /* Fluctuations for one eigenvector */
    evec = pczevec (pczf, iv);
    for (i = 0; i < pczf->natoms; i++) {
      j = 3 * i;
      fluc[i] = evec[j]*evec[j] + evec[j+1]*evec[j+1] + evec[j+2]*evec[j+2];
      fluc[i] = (float)sqrt(fluc[i]);
    }
    sfree (evec);
  } else {
    /* Fluctuations for all the eigenvectors */
    /* New method, using eigenvectors */
    /*
    for (i = 0; i < pczf->natoms; i++) {
      fluc[i] = 0.0;
    }
    for (k = 0 ; k < pczf->nvecs; k++) {
      evec = pczevec (pczf, k);
      for (i = 0; i < pczf->natoms; i++) {
        j = 3 * i;
        fluc[i] += (float)sqrt (evec[j]*evec[j] + evec[j+1]*evec[j+1] + evec[j+2]*evec[j+2]);
      }
      sfree (evec);
    }
    */
    
    /* Old method, using cartesian coordinates */
    memset (fluc, 0, pczf->natoms * sizeof (float));
    avg = pczavg (pczf);
    for (i = 0; i < pczf->nframes; i++) {
      snapshot = pczsnap (pczf, i);
      for (j = 0; j < pczf->natoms; j++) {
        x = snapshot[3*j]   - avg[3*j];
        y = snapshot[3*j+1] - avg[3*j+1];
        z = snapshot[3*j+2] - avg[3*j+2];
        fluc[j] += x*x + y*y + z*z;
      }
      sfree (snapshot);
    }
    sfree (avg);
    for (j = 0; j < pczf->natoms; j++) {
      fluc[j] = (float)sqrt (fluc[j]/(pczf->nframes));
    }
  }

  return fluc;
}

/* This method retro-projects a snapshot from the eigenvectors and
 * projections stored in the file. Is the method used to retrieve
 * the cartesian atomic data from the eigenvectors */
float *pczsnap (PCZFILE *pczf, int iframe) {
  float *snap, *scores, *evec;
  int i, ivec;

  if ((iframe < 0) || (iframe >= pczf->nframes))
    return NULL;

  snap = pczavg (pczf);
  if (snap == NULL)
    return NULL;

  scores = pczscores (pczf, iframe);
  if (scores == NULL)
    return NULL;

  for (ivec = 0; ivec < pczf->nvecs; ivec++) {
    evec = pczevec (pczf, ivec);
    if (evec == NULL)
      return NULL;

    for (i = 0; i < pczf->natoms * 3; i++) {
      snap[i] += evec[i] * scores[ivec];
    }
    sfree (evec);
  }
  sfree (scores);

  return snap;
}

/* This method transforms the atomic fluctuation data in the corresponding
 * B-factors, storing the result in the input array with the original data */
float *pczbfactors (float *fluc, int nAtoms) {
  int i;

  for (i = 0; i < nAtoms; i++)
    fluc[i] = (fluc[i] * fluc[i]) * ((8.0/3.0) * M_PI * M_PI);

  return fluc;
}

/* This method computes the Lindemann coefficient for the selected
 * atoms */
float pczlindemann (PCZFILE *pczf, int nAtoms, int *atomList) {
  float *fluc, *heavyfluc;
  int i;
  float lindemann;

  fluc = pczfluc (pczf, -1);
  lindemann = 0.0;
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
    heavyfluc = fluc;
  } else {
    heavyfluc = packValues (fluc, atomList, nAtoms);
    /*
    heavyfluc = (float *)smalloc (nAtoms * sizeof (float));
    for (i = 0; i < nAtoms; i++)
      heavyfluc[i] = fluc[atomList[i]];
    */
  }
  for (i = 0; i < nAtoms; i++)
    lindemann += heavyfluc[i];

  lindemann = (float)sqrt(lindemann / nAtoms) / 4.5;

  sfree (fluc);
  if (atomList != NULL)
    sfree (heavyfluc);

  return lindemann;
}

/* This method computes the foce constants for the selected
 * atoms */
float *pczforcecte (PCZFILE *pczf, float t) {
  float *evals, *ctes;
  int i;

  ctes = (float *)smalloc (pczf->nvecs * sizeof (float));
  evals = pczevals (pczf);

  for (i = 0; i < pczf->nvecs; i++) {
    if (evals[i] <= 0) {
      ctes[i] = 0.0;
    } else {
      ctes[i]= (float)(((R_CTE * t) / (F_CTE * evals[i])) * 1000.0); /* In cal/mol.A^2 */
    }
  }

  sfree (evals);

  return ctes;
}

/* This method returns the cartesian coordinates of a requested frame
 * with the requested eigenvector */
float *getFrameCoordinates (PCZFILE *pczf, int frameNumber, int nEvec) {
  float *snapshot, *scores, *evec;
  int j;

  if (nEvec == -1) {
    snapshot = pczsnap (pczf, frameNumber);
  } else {
    evec = pczevec (pczf, nEvec);
    snapshot = pczavg (pczf);
    scores = pczscores (pczf, frameNumber);
    for (j = 0; j < pczf->natoms * 3; j++) {
      snapshot[j] += evec[j] * scores[nEvec];
    }
    sfree (scores);
    sfree (evec);
  }

  return snapshot;
}
