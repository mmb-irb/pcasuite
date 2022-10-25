#ifndef VECTOR_H_
#define VECTOR_H_

typedef float Real;
typedef int Boolean;

void Vuv(Real *v);
void Vdiff(Real *dv, float *u, Real *v);
Real Vd2(Real *v1, Real *v2);
void Vvxv(Real *vxv, Real *v1, Real *v2);
void Vmv(Real *mv, Real *m, Real *v);
void Vmu(Real *mv, Real *m, float *v);
void Vmtv(Real *mv, Real *m, Real *v);
void Vtmtv(Real *mv, Real *m, float *v);
void Vmm(Real *mm, Real *m1, Real *m2);
void Vmmt(Real *mm, Real *m1, Real *m2);
void Vtmtv2(Real *mv, Real *m, float *v);

#endif /*VECTOR_H_*/
