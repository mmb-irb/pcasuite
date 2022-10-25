#ifndef PCAZIP_H_
#define PCAZIP_H_

#include <stdio.h>
#include "constants.h"

#define IDQ 90

struct clparameters_t {
  char *infile;
  char *outfile;
  char *pdbfile;
  char *maskfile;
  char *maskString;
  int natom;
  int nvecs;
  int iqual;
  int rms;
  int formatted;
  int verbose;
};
typedef struct clparameters_t CLParameters;

void calcproj (int natom3, int nframes, int nvecs, float *x, float *xt,
	       double *z, float *pn);
void cmrestore (double *cm, double *cmd, int n);
void cmsave (double *cm, double *cmd, int n);
void calccovar (float *x, float *xt, double *cm, int natom3, int nframes);
void pack (float *dest, float *source, int *mask, int length);
void packPdbAtom (PDBATOM *dest, PDBATOM *source, int *mask, int length);
void maskread (char *maskfile, int natom, int *nmask, int *mask);
int atomsInPDB (char *pdbfile);
void helptext (void);
char *copyString (char *str);
CLParameters *readCLParameters (int argc, char **argv);
void freeCLParameters (CLParameters *clp);
void checkCLP (CLParameters *clp);
double *computeEigenvalues (double *cm, int natom3);
double *computeEigenvectors (double *cm, int natom3, int nvecs);
int selectEigenvaluesByQuality (double *w, int length, float quality);
double getTotalVariance (double *w, int length);
double getVarianceForNEigenvalues (double *w, int length, int nEigenvalues);
double getQualityForNEigenvalues (double *w, int length, int nEigenvalues);

#ifdef __xlc__
#define dsyevr_ dsyevr
#endif

void dsyevr_ (char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA,
	double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W,
	double *Z, int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, int *IWORK,
	int *LIWORK, int *INFO);

#endif
