#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaussianrms.h"
#include "vector.h"
#include "constants.h"
#include "utils.h"

void gaussianfit (Real *xyz0, Real *xyz1, int natoms, int nframes) {
  int a, actFrame;
  float *wam; /* MXNA */
  float aux, *nextMoving;

  aux = 1.0 / natoms;
  /*nextMoving = xyz1;*/
  wam = (float *)smalloc (natoms * sizeof (float));
  for (actFrame = 0; actFrame < nframes; actFrame++) {
    nextMoving = xyz1 + (natoms * 3 * actFrame);
    /* Initializing Whole Weighting Vector */
    for (a = 0; a < natoms; a++) {
      wam[a] = aux;
    }

    /*
    if (verbose)
      printf ("Computing gaussian RMSd for frame %i\n", actFrame+1);
    */
    MHmatchWRms(xyz0, nextMoving, natoms, wam);
    /*nextMoving += natoms * 3;*/
  }
  sfree (wam);
}

void  MHmatchWRms(Real *xyz0, Real *xyz1, int natoms, Real *wam) {
  int     i, j, r;
  float   *rm;
  Real    v[3];
  int     a;
  float   *px0, *px1, *pw;
  Real    sumw;
  Real    d, rmsAt, oldRms, newRms, sumw_new, sumw_cull, rmsOrig;
  int     iter = 99;
  Real    trf[15];
  float   rms;
  float   *rmsA;
  float   EPS2 = 0.001;
  float   rmsAnt = 99.9;

  i = 1;
  rmsA = (float *)smalloc (natoms * sizeof (float));
  while (i <= iter) {
    /*printf ("Iter %d\n", i);*/
    for (j=0 ; j<15; j++)
      trf[j] = 0.0;
    KabschModified (natoms, xyz0, xyz1, wam, trf, TRUE, TRUE);
    px0 = xyz0;
    px1 = xyz1;
    pw = wam;
    rm = rmsA;
    sumw = sumw_new = sumw_cull = oldRms = newRms = 0;

    for (r = 0; r < natoms; r++) {
      Vtmtv(v, trf, px1);
      d = v[0] - *px0++;
      rmsAt = d * d;
      d = v[1] - *px0++;
      rmsAt += d * d;
      d = v[2] - *px0++;
      rmsAt += d * d;
      oldRms += rmsAt;
      *rm = rmsAt;
      newRms += *pw * rmsAt;
      sumw += *pw;
      /* New atom weight */
      /* *pw = 1 / (rmsAt + EPS2);*/
      *pw = exp (-rmsAt/2.0);

      sumw_new += *pw;
      pw++;
      px1 += 3;
    }
    rms =  sqrt (newRms / sumw);
    rmsOrig = sqrt (oldRms / natoms);

    pw = wam;
    for (a = 0; a < natoms; a++, pw++) {
      *pw = *pw/sumw_new;
    }

    if ((rmsAnt - rms) < EPS2 || i == iter) {
      /*
      if (verbose)
      	printf ("Kabsch End at iteration %d: %3.8f - %3.8f < %3.7f\n", i, rmsAnt, rms, EPS2);
      */
      break;
    }
    rmsAnt = rms;

    px1 = xyz1;
    for (a = 0; a < natoms; a++) {
      Vtmtv(v, trf, px1);
      *px1=v[0];px1++;
      *px1=v[1];px1++;
      *px1=v[2];px1++;
    }

    i++;
  }

	sfree (rmsA);
  /*
	if (verbose)
  	printf(" rmsFinal = %3.7f for %d atom pairs.\n\n", rms, natoms);
  */
}

/*-------- KabschModified: New Weighted RMS Superposition -------------*/
/*-------- Kabsch Modified 21-10-2005 ---------------------------------*/

/*------------------------------ NewTotalRMS = -----------------------------------------*/
/*------- = ( Sum(i) {RmsAt(i) * AtWeight(i)} ) / Sum(i) {AtWeight(i)}     -------------*/
/*-------     Where:                                                       -------------*/
/*-------       Sum(i)  =  Summatory (i = 0 ==> i = Number of Atoms)       -------------*/
/*-------       RmsAt  =  sqrt ((x-x1)^2 + (y-y1)^2 + (z-z1)^2)            -------------*/
/*-------       Atom Weight  =  1/(rmsAt*rmsAt)                            -------------*/
/*-------       By default, at First Iteration, All Atoms have Weight = 1  -------------*/
/*-------       Atom Weight Window of 5 (Optional)                         -------------*/

void KabschModified (int na, float xyz0[], float xyz1[], float w[], Real op[], Boolean cen0, Boolean cen1) {
/* int na;               number of atoms
 * float xyz0[], xyz1[]; reference (0) and moving (1) coordinate sets
 * float w[];            weight on each position
 * Real op[];            3x5 transformation operator (1->0)
 * Boolean cen0, cen1;   re-determine centr of ref, second str. if true
 */
  float *px0, *px1, *pw;
  Real  mxval, t;
  Real  R[3][3], RtR[3][3];
  Real  d[3], u[3]; /* eigenvalues */
  Real  e[3][3], v[9];  /* eigenvectors */
  Real  b[9];   /* = R * v */
  int   i, j, k, a;
  Real  aux;

  /* Remove translation */
  if (cen0) {
    Kcenter(na, xyz0, w, op);
    px0 = xyz0;
    for (a = 0; a < na; a++) {
      *px0++ -= op[0];
      *px0++ -= op[1];
      *px0++ -= op[2];
    }
  }
  if (cen1) {
    Kcenter(na, xyz1, w, op + 12);
    px1 = xyz1;
    for (a = 0; a < na; a++) {
      *px1++ -= op[12];
      *px1++ -= op[13];
      *px1++ -= op[14];
    }
  }

  /* Setup R and R'R matrices */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      R[i][j] = 0.0;
    }
  }
  px0 = xyz0;
  px1 = xyz1;
  pw = w;
  for (a = 0; a < na; a++, pw++) {
    for (i = 0; i < 3; i++) {
      aux = *pw * *(px0+i);
      for (j = 0; j < 3; j++) {
        R[i][j] += aux * *(px1+j);
      }
    }
    px0 += 3; px1 += 3;
  }
  mxval = 0.0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if ((t = fabs(R[i][j])) > mxval)
        mxval = t;
    }
  }
  if (mxval > 10.0)  {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        R[i][j] /= mxval;
      }
    }
  }
  Vmmt ((float *)RtR, (float *)R, (float *)R);  /* m' * m since 2D matrix is transpose of 1D */

  /* Get eigenvalues and eigenvectors of R'R */
  Kjacobi (RtR, d, e);

  /* and sort */
  k = 0;
  for (i = 1; i < 3; i++)
    if (d[i] > d[k])
      k = i;
  switch (k) {
    case 0:
      if (d[1] > d[2]) {j = 1; i = 2;}
      else {j = 2; i = 1;}
      break;
    case 1:
      if (d[0] > d[2]) {j = 0; i = 2;}
      else {j = 2; i = 0;}
    break;
    case 2:
      if (d[0] > d[1]) {j = 0; i = 1;}
      else {j = 1; i = 0;}
    break;
  }
  u[0] = d[k]; for (a = 0; a < 3; a++) v[a] = e[a][k];
  u[1] = d[j]; for (a = 0; a < 3; a++) v[a + 3] = e[a][j];
  u[2] = d[i]; Vvxv(v + 6, v, v + 3);

  /* Determine b = R * v */
  Vmtv(b, (float *)R, v);
  Vuv(b);
  Vmtv(b + 3, (float *)R, v + 3);
  Vuv(b + 3);
  Vvxv(b + 6, b, b + 3);

  /* The transformation matrix m = b*v' */
  Vmmt(op + 3, b, v);

  /* Remove translation removal */
  if (cen0) {
    px0 = xyz0;
    for (a = 0; a < na; a++) {
      *px0++ += op[0];
      *px0++ += op[1];
      *px0++ += op[2];
    }
  }
  if (cen1) {
    px1 = xyz1;
    for (a = 0; a < na; a++) {
      *px1++ += op[12];
      *px1++ += op[13];
      *px1++ += op[14];
    }
  }
}

/*---------- Calculate the coordinates of the center of a group of atoms ----*/
void Kcenter (int n, Real xyz[], float w[], Real xyzc[]) {
  int a;
  Real  *x, *y, *z, sw;

  x = xyzc;
  y = xyzc + 1;
  z = xyzc + 2;

  *x = 0.0;
  *y = 0.0;
  *z = 0.0;
  sw = 0.0;

  for (a = 0; a < n; a++) {
    *x += *w * *xyz++;
    *y += *w * *xyz++;
    *z += *w * *xyz++;
    sw += *w++;
  }

  *x /= sw;
  *y /= sw;
  *z /= sw;
}

/*--------------------------------------------------------------------
 * Kjacobi finds eigenvalues, d, and normalized eigenvectors, v,
 * of an NxN symmetric matrix a. N is currently set to 3.
 * This routine was initially written by Sue Lee in July 1991
 * by copying from the "Numerical Recipes in C".
 * On output, elements of a above the diagonal are destroyed.
 * BKL made small modifications in Jan. 1994, that consisted of introducing
 * the definition N and using Real type instead of the float type.
 * Real type can be typedef'ed to double or float.
 * -------------------------------------------------------------------- */
void Kjacobi(Real a[3][3], Real d[3], Real v[3][3]) {
  int j,iq,ip,i;
  Real  tresh,theta,tau,t,sm,s,h,g,c,b[3],z[3];
  int nrot;
  tresh = 1.0E-07;

  /* Initialize to the identity matrix */
  for (ip = 0; ip < 3; ip++) {
    for (iq = 0; iq < 3; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  /* Initialize b and d to the diagonal of a */
  /* This vector will accumulate terms of the form ta[pq] as in eq.11.1.14 */
  for (ip = 0; ip < 3; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  /* up to 50 Jacobi rotations */
  nrot = 0;
  for (i = 0; i < 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < 3-1; ip++) { /*sum off-diagonal elements*/
      for (iq = ip+1; iq < 3; iq++)
        sm += fabs(a[ip][iq]);
    }
    /* The normal return, which relies on quadratic convergence
     * to machine underflow */
    if (sm <= tresh)
      return;

    /* otherwise, cycle for each off-diagonal element */
    for (ip = 0; ip < 3-1; ip++) {
      for (iq = ip+1; iq < 3; iq++) {
        g = 100.0 * fabs (a[ip][iq]);
        /* After four sweeps, skip the rotation if the off-diagonal element is small */
        if ((i > 4) &&
          (((float)(fabs (d[ip]) + g)) == ((float)fabs (d[ip]))) &&
          (((float)(fabs (d[iq]) + g)) == ((float)fabs (d[iq]))))
          a[ip][iq] = 0.0;
        /* Normal return */


        if (fabs (a[ip][iq]) < tresh)
          return;
        /* Zero ip, iq element */
        h = d[iq] - d[ip];
        if (((float)(fabs (h) + g)) == ((float)fabs (h)))
          t = (a[ip][iq]) / h; /* t = 1/(2(theta))) */
        else {
          /* Equation 11.1.10 */
          theta = 0.5 * h / (a[ip][iq]);
          t = 1.0 / (fabs (theta) + sqrt (1.0 + (theta * theta)));
          if (theta < 0.0)
            t = -t;
        }
        c = 1.0 / sqrt (1 + t * t);
        s = t * c;
        tau = s / (1.0 + c);
        h = t * a[ip][iq];
        z[ip] -= h;
        z[iq] += h;
        d[ip] -= h;
        d[iq] += h;
        a[ip][iq] = 0.0;
        /* Case of rotations 1<j<p */
        for (j = 0; j < ip; j++) {
          /* Rotate(a,j,ip,j,iq) */
          g = a[j][ip];
          h = a[j][iq];
          a[j][ip] = g - s * (h + g * tau);
          a[j][iq] = h + s * (g - h * tau);
        }
        /* Case of rotations p<j<q */
        for (j = ip+1; j < iq; j++) {
          /* Rotate(a,ip,j,j,iq) */
          g = a[ip][j];
          h = a[j][iq];
          a[ip][j] = g - s * (h + g * tau);
          a[j][iq] = h + s * (g - h * tau);
        }
        /* Case of rotations q<j<3 */
        for (j = iq+1; j < 3; j++) {
          /* Rotate(a,ip,j,iq,j) */
          g = a[ip][j];
          h = a[iq][j];
          a[ip][j] = g - s * (h + g * tau);
          a[iq][j] = h + s * (g - h * tau);
        }
        for (j = 0; j < 3; j++) {
          /* Rotate(v,j,ip,j,iq) */
          g = v[j][ip];
          h = v[j][iq];
          v[j][ip] = g - s * (h + g * tau);
          v[j][iq] = h + s * (g - h * tau);
        }
        nrot++;
      }
    }
    for (ip = 0; ip < 3; ip++) {
      b[ip] += z[ip];
      /* Update d with the sum of ta[pq] and reinitialize z */
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  printf ("Too many iterations in routine JACOBI\n");
}
