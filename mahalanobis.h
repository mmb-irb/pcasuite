#ifndef MAHALANOBIS_H_
#define MAHALANOBIS_H_

#include <stdio.h>
#include "pcz_io.h"

void pczmahalanobis (PCZFILE *pczf, FILE *fout, int nvuse, char *xfile, int *atomList, int nAtoms);

#endif /*MAHALANOBIS_H_*/
