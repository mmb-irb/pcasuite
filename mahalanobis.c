#include <stdio.h>
#include <math.h>
#include "mahalanobis.h"
#include "pcz_io.h"
#include "traj_io.h"
#include "x_io.h"
#include "trajutils.h"
#include "utils.h"

/* This method computes the Mahalanobis distance for the reference data */
void pczmahalanobis (PCZFILE *pczf, FILE *fout, int nvuse, char *xfile, int *atomList, int nAtoms) {
  float *xavg, *evals, **evecs, *x, *proj, *aux, *xorig;
  float sm;
  int nvecs, nframes, i, j;
  TRAJFILE xf;

  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  }
  /* Reading of PCZ values */
  nvecs = (nvuse < pczf->nvecs) ? nvuse : pczf->nvecs;
  xavg = pczavg (pczf);
  evals = pczevals (pczf);
  if (atomList != NULL) {
    aux = packAtoms (xavg, atomList, nAtoms);
    sfree (xavg);
    xavg = aux;
  }
    
  evecs = (float **)smalloc (nvecs * sizeof (float *));
  for (i = 0; i < nvecs; i++) {
    evecs[i] = pczevec (pczf, i);
    evals[i] *= pczf->vtot/pczf->vexplained;
    if (atomList != NULL) {
      aux = packAtoms (evecs[i], atomList, nAtoms);
      sfree (evecs[i]);
      evecs[i] = aux;
    }
  }
  
  /* Reading of structures from trajectory */
  trajopen (xfile, pczf->natoms, &xf);
  nframes = xf.nframes;
  
  /* Distance computing */
  xorig = (float *)smalloc (pczf->natoms * 3 * sizeof (float));
  for (i = 1; i <= nframes; i++) {
    xsnap (&xf, i, xorig); /* Read frame */
    if (atomList != NULL) {
      x = packAtoms (xorig, atomList, nAtoms);
    } else {
      x = xorig;
    }
    trajfit (x, xavg, nAtoms, 1); /* Recenter the frame over the mean structure */
    for (j = 0; j < nAtoms * 3; j++) { /* Substract mean structure from recentered structure */
      x[j] -= xavg[j];
    }
    proj = matmul (evecs, x, nvecs, nAtoms * 3); /* Compute the projections */
    sm = 0.0;
    for (j = 0; j < nvecs; j++) { /* Compute the squared distances */
      sm += proj[j] * proj[j] / evals[j];
    }
    sm = sqrt (sm); /* Compute the distance */
    sfree (proj);
    if (atomList != NULL) {
      sfree (x);
    }
    
    fprintf (fout, "%8.3f\n", sm);
  }
  
  /* Resource freeing */
  fclose (xf.file);
  sfree (xorig);
  for (i = 0; i < nvecs; i++) {
    sfree (evecs[i]);
  }
  sfree (evecs);
  sfree (evals);
  sfree (xavg);
}
