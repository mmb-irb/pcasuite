#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include "pcz_io.h"
#include "constants.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define SIGN(x) ((x)>0.0?1.0:-1.0)

int fexist (char *filename);
char *copyString (char *str);
float minval (float *v, int nElems);
float maxval (float *v, int nElems);
char *trim (char *b);
void writeVector (float *vector, int nelems, FILE *fout);
void writeVectorD (double *vector, int nelems, FILE *fout);
void writeColumnVector (float *vector, int nelems, FILE *fout);
void writeCSV (float *vector, int nelems, FILE *fout);
void writeCSVSci (float *vector, int nelems, FILE *fout);
float *packAtoms (float *xyz, int *atomList, int nAtoms);
float *packValues (float *values, int *atomList, int nAtoms);
void computeMeanVariance (int *values, int nValues, float *mean, float *variance);
int am_big_endian (void);
float swap_endiannessf (float f);
unsigned int swap_endianness (unsigned int a);
char *rotateLeft(char *name, int length);
char *rotateRight(char *name, int length);
void *smalloc (size_t size);
void sfree (void *ptr);
void dumppdb (PCZFILE *pczf, FILE *fout, float *xyz, float *bfactors, int *atomList, int nAtoms);
float *matmul (float **matrix, float *vector, int rows, int columns);
long int getTime (void);
void pdbread (char *pdbfile, int natom, PDBATOM *pdba);
int addToVector (int *vector, int elem, int next);
int cmpintegers(const void *p1, const void *p2);

#endif
