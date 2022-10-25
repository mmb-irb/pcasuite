#ifndef GAUSSIANRMS_H_
#define GAUSSINARMS_H_

#include "vector.h"

typedef struct {
	Boolean	serial;
	int	bsn0;
	int bsn1;
	int	n;
	float	w;
} list_t;

void gaussianfit(Real *xyz0, Real *xyz1, int natoms, int nframes);
void MHmatchWRms(Real *xyz0, Real *xyz1, int natoms, Real *wam);
void KabschModified (int na, float xyz0[], float xyz1[], float w[], Real op[], Boolean cen0, Boolean cen1);
void Kcenter (int n, Real xyz[], float w[], Real xyzc[]);
void Kjacobi(Real a[3][3], Real d[3], Real v[3][3]);

#endif
