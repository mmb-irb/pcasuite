#ifndef PCZ_COMPUTE_H_
#define PCZ_COMPUTE_H_

#include "pcz_io.h"

float *pczfluc (PCZFILE *pczf, int iv);
float *pczsnap (PCZFILE *pczf, int iframe);
float *pczbfactors (float *fluc, int nAtoms);
float pczlindemann (PCZFILE *pczf, int nAtoms, int *atomMask);
float *pczforcecte (PCZFILE *pczf, float t);
float *getFrameCoordinates (PCZFILE *pczf, int frameNumber, int nEvec);

#endif
