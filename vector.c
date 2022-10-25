/*******************************************************************

* Program Kabsch Modified.
* New Functions:
- Weighted RMS (kabsch function modified).
- Amber Trajectories and Topologies Input.
- Amber Trajectories Output.
- Superposition of trajectories.
- Pdb Output with: Atom weights (Occupancy Columns).
Atom RMS (Temperature Factor Columns).
* Version Date: 21-10-2005
* Author: Adam Hospital Gasch (adam@mmb.pcb.ub.es).

********************************************************************/

/*---------------------------------*/
/*  Vector Manipulation Routines   */
/*---------------------------------*/

#include <math.h>
#include "vector.h"

#define	SMALL	1.0E-07

/*------- normalize a vector ------------------------------*/
void Vuv(Real *v) {
	Real	d;
	
	d = (Real)sqrt(*v * *v + *(v+1) * *(v+1) + *(v+2) * *(v+2));
	if (d < SMALL) {
		*v++ = 1.0;
		*v++ = 0.0;
		*v   = 0.0;
	}
	else {
		*v++ /= d;
		*v++ /= d;
		*v   /= d;
	}
}
/*------- dv = u - v, where u is of type float -------------*/
void Vdiff(Real *dv, float *u, Real *v) {
	*dv++ = *u++ - *v++;
	*dv++ = *u++ - *v++;
	*dv   = *u   - *v  ;
}
/*------- return |v1 - v2| * |v1 - v2| ---------------------*/
Real Vd2(Real *v1, Real *v2) {
	Real	d, d2;
	
	d = *v1++ - *v2++; d2  = d * d;
	d = *v1++ - *v2++; d2 += d * d;
	d = *v1   - *v2  ; d2 += d * d;
	return(d2);
}
/*----------vxv = v1 cross v2-------------------------------*/
void Vvxv(Real *vxv, Real *v1, Real *v2) {
	*vxv++ = *(v1+1) * *(v2+2) - *(v1+2) * *(v2+1);
	*vxv++ = *(v1+2) * *v2     - *v1     * *(v2+2);
	*vxv   = *v1     * *(v2+1) - *(v1+1) * *v2;
}
/*------- mv = m * v --------------------------------------*/
void Vmv(Real *mv, Real *m, Real *v) {
	*mv++ = *m     * *v + *(m+3) * *(v+1) + *(m+6) * *(v+2);
	*mv++ = *(m+1) * *v + *(m+4) * *(v+1) + *(m+7) * *(v+2);
	*mv   = *(m+2) * *v + *(m+5) * *(v+1) + *(m+8) * *(v+2);
}
/*------- like Vmv, except that vector v is of type float -----*/
void Vmu(Real *mv, Real *m, float *v) {
	*mv++ = *m     * *v + *(m+3) * *(v+1) + *(m+6) * *(v+2);
	*mv++ = *(m+1) * *v + *(m+4) * *(v+1) + *(m+7) * *(v+2);
	*mv   = *(m+2) * *v + *(m+5) * *(v+1) + *(m+8) * *(v+2);
}
/*------- mv = (m transpose) * v --------------------------------------*/
void Vmtv(Real *mv, Real *m, Real *v) {
	*mv++ = *m     * *v + *(m+1) * *(v+1) + *(m+2) * *(v+2);
	*mv++ = *(m+3) * *v + *(m+4) * *(v+1) + *(m+5) * *(v+2);
	*mv   = *(m+6) * *v + *(m+7) * *(v+1) + *(m+8) * *(v+2);
}
/*------- tmt is a 3x5 matrix of t0, m, and t1-------------
------- v' = tmt * v = t0 + m * (v - t1) ----------------
------- v is of type float ------------------------------*/
void Vtmtv(Real *mv, Real *m, float *v) {
	Real	*t, u[3];
	int	i;
	
	t = m + 12;
	for (i = 0; i < 3; i++) u[i] = *v++ - *t++;
	
	*mv++ = *m     + *(m+3) * *u + *(m+6) * *(u+1) + *(m+9) * *(u+2);
	*mv++ = *(m+1) + *(m+4) * *u + *(m+7) * *(u+1) + *(m+10) * *(u+2);
	*mv   = *(m+2) + *(m+5) * *u + *(m+8) * *(u+1) + *(m+11) * *(u+2);
}
/*------- mm = m1 * m2 --------------------------------------*/
void Vmm(Real *mm, Real *m1, Real *m2) {
	int	i;
	
	for (i = 0; i < 3; i++, m2 += 3) {
		*mm++ = *m1     * *m2 + *(m1+3) * *(m2+1) + *(m1+6) * *(m2+2);
		*mm++ = *(m1+1) * *m2 + *(m1+4) * *(m2+1) + *(m1+7) * *(m2+2);
		*mm++ = *(m1+2) * *m2 + *(m1+5) * *(m2+1) + *(m1+8) * *(m2+2);
	}
}
/*------- mm = m1 * (m2 transpose) ---------------------------*/
void Vmmt(Real *mm, Real *m1, Real *m2) {
	int	i;

	for (i = 0; i < 3; i++, m2++) {
		*mm++ = *m1     * *m2 + *(m1+3) * *(m2+3) + *(m1+6) * *(m2+6);
		*mm++ = *(m1+1) * *m2 + *(m1+4) * *(m2+3) + *(m1+7) * *(m2+6);
		*mm++ = *(m1+2) * *m2 + *(m1+5) * *(m2+3) + *(m1+8) * *(m2+6);
	}
}
/*------- tmt is a 3x5 matrix of t0, m, and t1-------------
------- v is of type float ------------------------------
------- Difference between Vtmtv and Vtmtv2 -------------
------- is the translation factor, here there's not -----
------- translation factor. -----------------------------*/
void Vtmtv2(Real *mv, Real *m, float *v) {
	*mv++ = *m     + *(m+3) * *v + *(m+6) * *(v+1) + *(m+9) * *(v+2);
	*mv++ = *(m+1) + *(m+4) * *v + *(m+7) * *(v+1) + *(m+10) * *(v+2);
	*mv   = *(m+2) + *(m+5) * *v + *(m+8) * *(v+1) + *(m+11) * *(v+2);
}
/*------- end vector routines --------------*/
