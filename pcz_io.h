#ifndef PCZ_IO_H_
#define PCZ_IO_H_

#include <stdio.h>
#include "constants.h"

/* This definitions serve the purpose of calculate the space needed by every
 * data scratch depending on the file version number. They could be methods,
 * but this approach is faster. Either case, care must be taken when using
 * this macros in order to avoid evaluating pczf aprameter multiple times. */
#define HDR_SIZE(pczf) ((pczf)->pczVersion==4?(116+((pczf)->haveAtomNames?(16*((pczf)->natoms)):0)):(100))
#define MEAN_SIZE(pczf) (4*3*((pczf)->natoms))
#define EVEC_SIZE(pczf) (4*3*((pczf)->natoms))
#define EVAL_SIZE(pczf) 4
#define PROJ_SIZE(pczf) (4*((pczf)->nframes))
#define EIGEN_DATA_SIZE(pczf) (EVEC_SIZE(pczf)+EVAL_SIZE(pczf)+PROJ_SIZE(pczf))
#define PADDING(pczf) ((pczf)->pczVersion==2?(EIGEN_DATA_SIZE(pczf)-(HDR_SIZE(pczf)+MEAN_SIZE(pczf))):0)

struct pczfile_t {
  char filename[81];
  char title[81];
  char version[4];
  int pczVersion;
  FILE *file;
  int status;
  int natoms;
  int nvecs;
  int nframes;
  float vtot;
  float vexplained;
  float quality;
  int dimensionality;
  int rms;
  int haveAtomNames;
	PDBATOM *atomNames;
};
typedef struct pczfile_t PCZFILE;

PCZFILE *pczopen (char *fname);
void pczclose (PCZFILE *pczf);
float *pczevals (PCZFILE *pczf);
float *pczevec (PCZFILE *pczf, int nEvec);
float *pczproj (PCZFILE *pczf, int ip);
float *pczavg (PCZFILE *pczf);
float *pczfluc (PCZFILE *pczf, int iv);
float *pczscores (PCZFILE *pczf, int ir);
float *pczbfactors (float *fluc, int nAtoms);
float pczlindemann (PCZFILE *pczf, int nAtoms, int *atomMask);
float *pczforcecte (PCZFILE *pczf, float t);

void pczwrite (char *fName, PCZFILE *pczf, float *mean, float *evals, float **evecs, float **proj);

#endif
