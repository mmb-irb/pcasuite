#ifndef MPI_ROUTINES_H_
#define MPI_ROUTINES_H_

long usecs (void);
float *computeSum (float *array, int nrows, int ncols);
float *computeMean (float *array, int nrows, int ncols, int nelems);
void meanCovar (float *A, int filas, int columnas, float *mean, double **outCovar);
int min (int a, int b);
void computeEigen (double *inMatrix, int N, double *quality, double **evals, double **evecs, int *howManyVectors, float *vtot, float *vsum, int *dimensionality);
void calcproj_mpi (int inColumns, int inNframes, int inNvecs, float *inXyz, float *inMean, double *inEvecs, float *outProj);
void trajfit_mpi (float *x, float *ref, int natoms, int nframes);
void gaussianfit_mpi (float *inXyz0, float *inXyz1, int inNatoms, int inNframes);
void getDataDistribution (int n, int howMany, int **outFilasEnCPU, int **outDesplParaCPU);

#endif /*MPI_ROUTINES_H_*/
