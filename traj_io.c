#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "traj_io.h"
#include "pcazip.h"
#include "utils.h"
#include "x_io.h"
#include "binpos_io.h"
#include "netcdf_io.h"

/* The next three arrays must be coherent: the n-th position of any array
 * must contain the adequate function of a format module, i.e. functions
 * in the n-th position must be able to deal with binpos format */
/***************************************************************************/
/* Array with the pointers to the functions that open files in a particular
 * format */
trjopen_t trjopen[] = {
  &binposopen,
  &netcdfopen,
  &xopen,
  NULL};
/* Array with the pointers to the functions that read files in a particular
 * format */
trjsnap_t trjsnap[] = {
  &binpossnap,
  &netcdfsnap,
  &xsnap,
  NULL};
/* Array with the pointers to the functions that check if a file can be read
 * with the corresponding functions */
trjFormatOK_t trjFormatOK[] = {
  &binposFormatOK,
  &netcdfFormatOK,
  &xFormatOK,
  NULL};
/* Array with the pointers to the functions that closes files in a particular
 * format */
trjopen_t trjclose[] = {
  &binposclose,
  &netcdfclose,
  &xclose,
  NULL};
/***************************************************************************/

/* This method search for a proper module to read the data */
void trajopen (char *fname, int natoms, TRAJFILE *xf) {
  int i, accepted;
  
  xf->filename = copyString (fname);
  xf->natoms = natoms;
  for (i = 0, accepted = FALSE; trjFormatOK[i] != NULL && accepted == FALSE; i++) {
    accepted = (*trjFormatOK[i]) (xf);
    if (accepted) {
      xf->trjopen  = trjopen[i];
      xf->trjsnap  = trjsnap[i];
      xf->trjclose = trjclose[i];
    }
  }
  if (accepted) {
    (*xf->trjopen) (xf);
  } else {
    fprintf (stderr, "error: Unrecognized input file format\n");
    exit (3);
  }
}

/* This method searchs for a given frame (is) into the file and returns
 * the coordinates stored in that frame in the supplied buffer */
int trajsnap (TRAJFILE *xf, int is, float *buffer) {
  return (*xf->trjsnap)(xf, is, buffer);
}

void trajclose (TRAJFILE *xf) {
  (*xf->trjclose)(xf);
  xf->file = NULL;
  xf->fd = 0;
}

void trajfree (TRAJFILE *xf) {
  if (xf->filename != NULL)
    sfree (xf->filename);
  if (xf->title != NULL)
    sfree (xf->title);
}

/* This method returns a boolean indicating if the trajectory has box */
int trajbox (TRAJFILE *xf) {
  return (xf->hasbox == TRUE);
}

/* This method returns the number of frames that the trajectory has */
int trajframes (TRAJFILE *xf) {
  return xf->nframes;
}

/* This method returns an array with all the coordinates and frames loaded */
float **trajmemread (TRAJFILE *xf) {
  float **x;
  int i;

  x = (float **)malloc (xf->nframes * sizeof (float *));
  for (i = 0; i < xf->nframes; i++) {
    x[i] = (float *)malloc (xf->natoms * 3 * sizeof (float));
    (*xf->trjsnap) (xf, i + 1, x[i]);
  }

  return x;
}
