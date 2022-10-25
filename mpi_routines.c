#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include "mpi_routines.h"
#include "trajutils.h"
#include "gaussianrms.h"
#include "scalapack.h"
#include "utils.h"

long usecs (void) {
  struct timeval t;

  gettimeofday(&t,NULL);

  return t.tv_sec*1000000+t.tv_usec;
}

float *computeSum (float *array, int nrows, int ncols) {
  float *sum;
  int i, j;

  sum = (float *)malloc (ncols * sizeof (float));
  bzero (sum, ncols * sizeof (float));
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      /*sum[j] += (double)(array[i*ncols+j]);*/
      sum[j] += array[j*nrows+i];
    }
  }

  return sum;
}

float *computeMean (float *array, int nrows, int ncols, int nelems) {
  float *mean;
  int i;

  mean = computeSum (array, nrows, ncols);
  for (i = 0; i < ncols; i++) {
    mean[i] /= (float)nelems;
  }

  return mean;
}

void meanCovar (float *A, int filas, int columnas, float *mean, double **outCovar) {
  int i, j, howMany, whoAmI, var1, var2, fila, columna;
  int *filasEnCPU, *desplParaCPU, filasPorCPU, filasPendientes, despl;
  double val1, val2;
  double *pCovar, *covar, **result;
  float *array, *sum;
  MPI_Datatype RowType;
  MPI_Aint addressOut, addressIn;

  MPI_Comm_size (MPI_COMM_WORLD, &howMany);
  MPI_Comm_rank (MPI_COMM_WORLD, &whoAmI);
  
  /* Distribute scalar input parameters to other processors */
  MPI_Bcast (&filas,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&columnas, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Prepare an MPI type for each row */
  MPI_Type_contiguous (columnas, MPI_FLOAT, &RowType);
  MPI_Type_commit (&RowType);
  
  /* Generate distribution of rows */
  getDataDistribution (filas, howMany, &filasEnCPU, &desplParaCPU);

  /* Prepare space for MPI transfers */
  array = (float *)smalloc ((filasEnCPU[whoAmI])*columnas*sizeof (float));

  /* Distribute data to processors */
  MPI_Address (&A[0], &addressOut);
  MPI_Address (&array[0], &addressIn);
  MPI_Scatterv ((void *)addressOut, filasEnCPU, desplParaCPU, RowType, (void *)addressIn, filasEnCPU[whoAmI], RowType, 0, MPI_COMM_WORLD);
  
  /* Free of used data */
  free (desplParaCPU);

  /* The mean has been precomputed, so distribute the input mean */
  if (whoAmI != 0) {
    /* Create the mean buffer to receive the broadcast */
    mean = (float *)smalloc (columnas * sizeof (float));
  }

  /* Redistribute the mean to all the processors */
  MPI_Address (&mean[0], &addressIn);
  MPI_Bcast ((void *)addressIn, 1, RowType, 0, MPI_COMM_WORLD);

  /* Now we compute the partial covariance matrices */
  pCovar = (double *)smalloc (columnas * columnas * sizeof (double));
  bzero (pCovar, columnas * columnas * sizeof (double));
  for (fila = 0; fila < filasEnCPU[whoAmI]; fila++) {
    for (var1 = 0; var1 < columnas; var1++) {
      val1  = (double)array[fila*columnas+var1];
      val1 -= (double)mean[var1];
      for (var2 = var1; var2 < columnas; var2++) {
        val2  = (double)array[fila*columnas+var2];
        val2 -= (double)mean[var2];
        pCovar[var1*columnas+var2] += val1 * val2;
      }
    }
  }

  /* Finally we send the matrices back for summation and final processing */
  if (whoAmI == 0) {
    covar = (double *)smalloc (columnas*columnas*sizeof (double));
  }
  MPI_Address (&pCovar[0], &addressOut);
  MPI_Address (&covar[0], &addressIn);
  MPI_Reduce ((void *)addressOut, (void *)addressIn, columnas * columnas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (whoAmI == 0) {
    for (fila = 0; fila < columnas; fila ++) {
      for (columna = fila; columna < columnas; columna ++) {
        covar[fila*columnas+columna] /= (double)filas;
        covar[columna*columnas+fila] = covar[fila*columnas+columna];
      }
    }
  }


  free (array);
  free (pCovar);
  free (filasEnCPU);
  MPI_Type_free (&RowType);

  if (whoAmI==0) {
    /* Result return in master node */
    *outCovar = covar;
  } else {
    /* Resource freeing in slave nodes */
    free (mean);
  }
}

int min (int a, int b) {
  return a<b?a:b;
}

void computeEigen (double *inputMatrix, int inN, double *inQuality, double **evals, double **evecs, int *inHowManyVectors, float *vtot, float *vsum, int *dimensionality) {
  char jobz, uplo;
  int ZERO = 0, ONE = 1;
  int iAm, nProcs, nb, nprow, npcol, myRow, myCol, myRows, myCols, ctxt, info;
  int descA[9], descZ[9], lwork;
  int i, j, iu, il;
  double *lA, *W, *Z, *lZ, *work;
  int N, howManyVectors;
  double quality, *inMatrix;

  /* Blocking factor for matrix distribution */
  nb = 64;

  /* BLACS initialization */
  blacs_pinfo (&iAm, &nProcs);

  /* Compute an optimum number of processors */
  nprow = (int)sqrt ((double)nProcs);
  npcol = nProcs / nprow;
  if ((nprow*npcol < nProcs) && (nProcs < 9)) {
    nprow = 1;
    npcol = nProcs;
  }

  blacs_get (NULL, &ZERO, &ctxt);
  blacs_gridinit (&ctxt, "R", &nprow, &npcol);
  blacs_gridinfo (&ctxt, &nprow, &npcol, &myRow, &myCol);
  
  if (myRow == -1 && myCol == -1) {
    /* I am out of the grid. No computation for me. */
    return;
  }

  /* Distribute the matrix dimension */
  if (myRow == 0 && myCol == 0) {
    /* We are the master node. Send the data */
    N = inN;
    quality = *inQuality;
    howManyVectors = *inHowManyVectors;
    inMatrix = inputMatrix;
    igebs2d (&ctxt, "A", " ", &ONE, &ONE, &N, &ONE);
    dgebs2d (&ctxt, "A", " ", &ONE, &ONE, &quality, &ONE);
    igebs2d (&ctxt, "A", " ", &ONE, &ONE, &howManyVectors, &ONE);
  } else {
    /* We are the slave nodes. Receive the data */
    igebr2d (&ctxt, "A", " ", &ONE, &ONE, &N, &ONE, &ZERO, &ZERO);
    dgebr2d (&ctxt, "A", " ", &ONE, &ONE, &quality, &ONE, &ZERO, &ZERO);
    igebr2d (&ctxt, "A", " ", &ONE, &ONE, &howManyVectors, &ONE, &ZERO, &ZERO);
  }

  /* How many rows and cols the current processor will have */
  myRows = numroc (&N, &nb, &myRow, &ZERO, &nprow);
  myCols = numroc (&N, &nb, &myCol, &ZERO, &npcol);

  /* Prepare the matrices, initialize descriptors and allocate memory */
  /* Local original matrix */
  lA = (double *)calloc (myRows*myCols, sizeof(double));
  /* Eigenvalues vector */
  W  = (double *)calloc (N, sizeof (double));
  /* Local eigenvectors matrix */
  lZ = (double *)calloc (myRows*myCols, sizeof (double));
  /* Descriptors */
  descinit (descA, &N, &N, &nb, &nb, &ZERO, &ZERO, &ctxt, &myRows, &info);
  descinit (descZ, &N, &N, &nb, &nb, &ZERO, &ZERO, &ctxt, &myRows, &info);

  /* Global matrix distribution to all processors. The matrix must be upper diagonal */
  if (myRow == 0 && myCol == 0) {
    /* We are the master node. Send the data */
    dtrbs2d (&ctxt, "A", " ", "U", "L", &N, &N, inMatrix, &N);
  } else {
    /* We are the slave nodes. Receive the data */
    inMatrix = (double *)calloc (N * N, sizeof (double));
    dtrbr2d (&ctxt, "A", " ", "U", "L", &N, &N, inMatrix, &N, &ZERO, &ZERO);
  }

  /* Local matrix filling with global data */
  for (i = 1; i <= N; i++) {
    for (j = i; j <= N; j++) {
      pdelset (lA, &i, &j, descA, &inMatrix[((j-1) * N) + (i-1)]);
    }
  }
  /* Diagonalization, the hard work */
  /* Workspace query */
  jobz = 'N';
  uplo = 'U';
  work = (double *)calloc (2, sizeof (double));
  lwork = -1;
  pdsyev (&jobz, &uplo, &N, lA, &ONE, &ONE, descA, W, lZ, &ONE, &ONE, descZ, work, &lwork, &info);

  /* Eigenvalues computation */
  lwork = (int)work[0];
  free (work);
  work = (double *)calloc (lwork, sizeof (double));
  pdsyev (&jobz, &uplo, &N, lA, &ONE, &ONE, descA, W, lZ, &ONE, &ONE, descZ, work, &lwork, &info);

  /* Freeing work areas */
  free (work);

  /* Compute which eigenvectors need to be calculed for the desired quality */
  *vtot = 0.0;
  for (i = 0; i < N; i++) {
    *vtot += W[i];
  }

  iu = N;
  *vsum = 0.0;
  if (quality <= 0.0 && howManyVectors <= 0) {
    quality = (double)90;
  }
  if (quality > 0.0) {
    for (i = N-1; i >= 0; i--) {
      *vsum += W[i];
      if ((*vsum / *vtot) > (quality/100.0)) {
        il = i + 1;
        break;
      }
    }
  } else {
    il = N - (howManyVectors) + 1;
    for (i = il - 1; i < N; i++) {
      *vsum += W[i];
    }
  }
  quality = *vsum / *vtot * 100.0;
  howManyVectors = iu - il + 1;


  /* Return value assignment for eigenvectors after dimensionality calculus */
  if (myRow == 0 && myCol == 0) {
    *dimensionality = 0;
    while ((N - (*dimensionality) - 1) >= 0 && W[N - (*dimensionality) - 1] >= 1.0) {
      /*fprintf (stderr, "Big dimensionality: %4i->%f\n", N - (*dimensionality) - 1, W[N - (*dimensionality) - 1]);*/
      (*dimensionality)++;
    }
    
    *evals = W;
  } else {
    free (W);
  }

  /* Eigenvalues computation destroyed the diagonal elements and
   * the upper part of the matrix. Refill the local arrays from the original
   * matrix for the eigenvectors computation */
  for (i = 1; i <= N; i++) {
    for (j = i; j <= N; j++) {
      pdelset (lA, &i, &j, descA, &inMatrix[((j-1) * N) + (i-1)]);
    }
  }

  /* Diagonalization, second part: Eigenvector computation */
  /* Workspace query */
  jobz  = 'V';
  char range = 'I';
  uplo  = 'U';
  double abstol = (double)0;
  W  = (double *)calloc (N, sizeof (double));
  int m; /* Total number of evals found */
  int nz; /* Total number of evecs computed */
  double orfac = (double)0;
  work = (double *)calloc (3, sizeof (double));
  lwork = -1;
  int *iwork = (int *)calloc (1, sizeof (int));
  int liwork = -1;
  int *ifail = (int *)calloc (N, sizeof (int));
  int *iclustr = (int *)calloc (2 * nprow * npcol, sizeof (int));
  double *gap = (double *)calloc (nprow * npcol, sizeof (double));
  double vl = 0.0, vu = 0.0;
  
  pdsyevx (&jobz, &range, &uplo, &N, lA, &ONE, &ONE, descA, &vl, &vu, &il, &iu, &abstol, &m, &nz, W, &orfac, lZ, &ONE, &ONE, descZ, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);
  
  /* Eigenvectors computation */
  lwork = (int)work[0];
  free (work);
  work = (double *)calloc (lwork, sizeof (double));
  liwork = iwork[0];
  free (iwork);
  iwork = (int *)calloc (liwork, sizeof (int));
  pdsyevx (&jobz, &range, &uplo, &N, lA, &ONE, &ONE, descA, &vl, &vu, &il, &iu, &abstol, &m, &nz, W, &orfac, lZ, &ONE, &ONE, descZ, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);

  /* Freeing work areas  and other not used vectors */
  free (work);
  free (iwork);
  free (W);
  free (ifail);
  free (iclustr);
  free (gap);

  /* Retrieving the distributed results back to a global matrix */
  /* Global eigenvectors matrix */
  Z  = (double *)calloc (N * (iu - il + 1), sizeof (double));
  for (i = 1; i <= N; i++) {
    for (j = 1; j <= iu - il + 1; j++) {
      pdelget ("A", " ", &Z[(j-1)*N+(i-1)], lZ, &i, &j, descZ);
    }
  }

  /* Return values */
  *evecs = Z;
  if (myRow == 0 && myCol == 0) {
    *inQuality = quality;
    *inHowManyVectors = howManyVectors;
  }



  /* Resource freeing */
  if (myRow != 0 || myCol != 0) {
    free (inMatrix);
    free (Z);
  }
  free (lA);
  free (lZ);
  blacs_gridexit (&ctxt);

  /* This routine is intended to be called from an MPI program which
   * uses MPI for other purposes, so we don't explicitly close the
   * MPI communications, leaving the job to the caller */
  blacs_exit (&ONE);
}

/* Method for projection calculus
 * Parameters:
 * columns -> (Input)  Number of atoms *3
 * nframes -> (Input)  Number of frames in trajectory
 * nvecs   -> (Input)  Number of eigenvectors
 * xyz     -> (Input)  Trajectory
 * mean    -> (Input)  Mean structure
 * evecs   -> (Input)  Eigenvectors
 * proj    -> (Output) Projection
 */
void calcproj_mpi (int inColumns, int inNframes, int inNvecs, float *inXyz,
  float *inMean, double *inEvecs, float *outProj) {
  int i, j, k, l, whoAmI, howMany;
  int columns, nframes, nvecs;
  int *filasEnCPU, *desplParaCPU, filasPorCPU, filasPendientes, despl;
  float *xyz, *mean, *proj, *lxyz, *lproj;
  double *evecs;
  MPI_Datatype RowType;
  MPI_Aint addressOut, addressIn;
  
  /* Get MPI information */
  MPI_Comm_size (MPI_COMM_WORLD, &howMany);
  MPI_Comm_rank (MPI_COMM_WORLD, &whoAmI);

  /* Initialize master variables */
  if (whoAmI == 0) {
    columns = inColumns;
    nframes = inNframes;
    nvecs   = inNvecs;
    xyz     = inXyz;
    mean    = inMean;
    evecs   = inEvecs;
    proj    = outProj;
  }
  
  /* Distribute scalar input parameters to all processors */
  MPI_Bcast (&columns, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&nframes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&nvecs,   1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Distribute vector data  to all processors */
  getDataDistribution (nframes, howMany, &filasEnCPU, &desplParaCPU);
  /* Get space for input arrays */
  if (whoAmI != 0) {
    mean  = (float *) smalloc (columns * sizeof (float));
    evecs = (double *)smalloc (columns * nvecs * sizeof (double));
  }
  /* Distribute whole-size arrays */
  MPI_Bcast (mean,  columns,         MPI_FLOAT,  0, MPI_COMM_WORLD);
  MPI_Bcast (evecs, columns * nvecs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  /* Get space for fragmented and local arrays */
  lxyz  = (float *)smalloc (filasEnCPU[whoAmI] * columns * sizeof (float));
  lproj = (float *)smalloc (filasEnCPU[whoAmI] * nvecs * sizeof (float));

  /* Prepare an MPI type for each row */
  MPI_Type_contiguous (columns, MPI_FLOAT, &RowType);
  MPI_Type_commit (&RowType);

  /* Distribute fragmented array */
  MPI_Address (&xyz[0], &addressOut);
  MPI_Address (&lxyz[0], &addressIn);
  MPI_Scatterv ((void *)addressOut, filasEnCPU, desplParaCPU, RowType, (void *)addressIn, filasEnCPU[whoAmI], RowType, 0, MPI_COMM_WORLD);

  /* Compute the projections for the data received */
  for (k = 0; k < filasEnCPU[whoAmI]; k++) {
    /* First substract average structure */
    for (i = 0; i < columns; i++) {
      lxyz[k * columns + i] -= mean[i];
    }
  
    /* Loop over eigenvectors */
    for (j = 0; j < nvecs; j++) {
      lproj[j * filasEnCPU[whoAmI] + k] = 0.0;
      for (l = 0; l < columns; l++) {
        lproj[j * filasEnCPU[whoAmI] + k] +=
          evecs[j * columns + l] * lxyz[k * columns + l];
      }
    }
  }

  /* Return the projections to master node */
  for (i=0; i<nvecs; i++) {
    MPI_Address (&lproj[i*filasEnCPU[whoAmI]], &addressOut);
    MPI_Address (&proj[i*nframes], &addressIn);
    MPI_Gatherv ((void *)addressOut, filasEnCPU[whoAmI], MPI_FLOAT, (void *)addressIn, filasEnCPU, desplParaCPU, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }
  
  /* Resource freeing */
  MPI_Type_free (&RowType);
  sfree (filasEnCPU);
  sfree (desplParaCPU);
  sfree (lproj);
  sfree (lxyz);
  if (whoAmI != 0) {
    sfree (mean);
    sfree (evecs);
  }
}

/* This method receives a trajectory (x), a reference frame (ref) and
 * modifies the trajectory to the RMSd best fit with the reference frame
 * using the Kabsch algorithm (standard version) */
void trajfit_mpi (float *inX, float *inRef, int inNatoms, int inNframes) {
  int findmove, i, j;
  float xt, yt, zt, rmsd;
  float r[9], v[3];
  float *x, *ref;
  int natoms, nframes;
  int whoAmI, howMany, *filasEnCPU, *desplParaCPU;
  MPI_Datatype RowType;
  MPI_Aint addressOut, addressIn;

  /* Get MPI information */
  MPI_Comm_size (MPI_COMM_WORLD, &howMany);
  MPI_Comm_rank (MPI_COMM_WORLD, &whoAmI);

  /* Initialize master variables */
  if (whoAmI == 0) {
    natoms  = inNatoms;
    nframes = inNframes;
    ref     = inRef;
  }
  
  /* Distribute scalar input parameters to all processors */
  MPI_Bcast (&natoms,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&nframes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Distribute vector data  to all processors */
  getDataDistribution (nframes, howMany, &filasEnCPU, &desplParaCPU);

  /* Prepare an MPI type for each row */
  MPI_Type_contiguous (natoms * 3, MPI_FLOAT, &RowType);
  MPI_Type_commit (&RowType);

  /* Prepare space for MPI transfers */
  x = (float *)smalloc ((filasEnCPU[whoAmI]) * natoms * 3 * sizeof (float));
  if (whoAmI != 0) {
    ref = (float *)smalloc (natoms * 3 * sizeof (float));
  }

  /* Distribute data to processors */
  /* Reference structure distribution to all */
  MPI_Address (&ref[0], &addressIn);
  MPI_Bcast ((void *)addressIn, 1, RowType, 0, MPI_COMM_WORLD);
  /* Trajectory structures distribution to each processor */
  MPI_Address (&inX[0], &addressOut);
  MPI_Address (&x[0], &addressIn);
  MPI_Scatterv ((void *)addressOut, filasEnCPU, desplParaCPU, RowType, (void *)addressIn, filasEnCPU[whoAmI], RowType, 0, MPI_COMM_WORLD);
  
  /* RMSd fitting */
  findmove = TRUE;
  
  for (i = 0; i < filasEnCPU[whoAmI]; i++) {
    matfit(&natoms, ref, &x[i*natoms*3], r, v, &rmsd, &findmove);
    for (j = 0; j < natoms; j++) {
      xt = r[0] * x[0+((j+(i*natoms))*3)] +
           r[3] * x[1+((j+(i*natoms))*3)] +
           r[6] * x[2+((j+(i*natoms))*3)] + v[0];
      yt = r[1] * x[0+((j+(i*natoms))*3)] +
           r[4] * x[1+((j+(i*natoms))*3)] +
           r[7] * x[2+((j+(i*natoms))*3)] + v[1];
      zt = r[2] * x[0+((j+(i*natoms))*3)] +
           r[5] * x[1+((j+(i*natoms))*3)] +
           r[8] * x[2+((j+(i*natoms))*3)] + v[2];
      x[0+((j+(i*natoms))*3)] = xt;
      x[1+((j+(i*natoms))*3)] = yt;
      x[2+((j+(i*natoms))*3)] = zt;
    }
  }

  /* Data gathering */
  MPI_Address (&x[0], &addressOut);
  MPI_Address (&inX[0], &addressIn);
  MPI_Gatherv ((void *)addressOut, filasEnCPU[whoAmI], RowType, (void *)addressIn, filasEnCPU, desplParaCPU, RowType, 0, MPI_COMM_WORLD);
  
  /* Resource freeing */
  sfree (filasEnCPU);
  sfree (desplParaCPU);
  sfree (x);
  if (whoAmI != 0) {
    sfree (ref);
  }
  
}

void gaussianfit_mpi (float *inXyz0, float *inXyz1, int inNatoms, int inNframes) {
  int a, actFrame;
  float *wam; /* MXNA */
  float aux, *nextMoving;
  int natoms, nframes;
  float *xyz0, *xyz1;
  int whoAmI, howMany, *filasEnCPU, *desplParaCPU;
  MPI_Datatype RowType;
  MPI_Aint addressOut, addressIn;

  /* Get MPI information */
  MPI_Comm_size (MPI_COMM_WORLD, &howMany);
  MPI_Comm_rank (MPI_COMM_WORLD, &whoAmI);

  /* Initialize master variables */
  if (whoAmI == 0) {
    natoms  = inNatoms;
    nframes = inNframes;
    xyz0    = inXyz0;
  }
  
  /* Distribute scalar input parameters to all processors */
  MPI_Bcast (&natoms,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&nframes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Distribute vector data  to all processors */
  getDataDistribution (nframes, howMany, &filasEnCPU, &desplParaCPU);

  /* Prepare an MPI type for each row */
  MPI_Type_contiguous (natoms * 3, MPI_FLOAT, &RowType);
  MPI_Type_commit (&RowType);

  /* Prepare space for MPI transfers */
  xyz1 = (float *)smalloc ((filasEnCPU[whoAmI]) * natoms * 3 * sizeof (float));
  if (whoAmI != 0) {
    xyz0 = (float *)smalloc (natoms * 3 * sizeof (float));
  }

  /* Distribute data to processors */
  /* Reference structure distribution to all */
  MPI_Address (&xyz0[0], &addressIn);
  MPI_Bcast ((void *)addressIn, 1, RowType, 0, MPI_COMM_WORLD);
  /* Trajectory structures distribution to each processor */
  MPI_Address (&inXyz1[0], &addressOut);
  MPI_Address (&xyz1[0], &addressIn);
  MPI_Scatterv ((void *)addressOut, filasEnCPU, desplParaCPU, RowType, (void *)addressIn, filasEnCPU[whoAmI], RowType, 0, MPI_COMM_WORLD);
  
  /* RMSd fitting */
  aux = 1.0 / natoms;
  wam = (float *)smalloc (natoms * sizeof (float));
  for (actFrame = 0; actFrame < filasEnCPU[whoAmI]; actFrame++) {
    nextMoving = xyz1 + (natoms * 3 * actFrame);
    /* Initializing Whole Weighting Vector */
    for (a = 0; a < natoms; a++) {
      wam[a] = aux;
    }

    MHmatchWRms(xyz0, nextMoving, natoms, wam);
  }
  sfree (wam);

  /* Data gathering */
  MPI_Address (&xyz1[0], &addressOut);
  MPI_Address (&inXyz1[0], &addressIn);
  MPI_Gatherv ((void *)addressOut, filasEnCPU[whoAmI], RowType, (void *)addressIn, filasEnCPU, desplParaCPU, RowType, 0, MPI_COMM_WORLD);

  /* Resource freeing */
  sfree (filasEnCPU);
  sfree (desplParaCPU);
  sfree (xyz1);
  if (whoAmI != 0) {
    sfree (xyz0);
  }

}


/* This method return two arrays for use in MPI_Gather/Scatter kind of
 * functions. Given the munber of elements (n) to distribute over a number
 * of processors (howMany), it returns an array with the number of elements
 * for each processor and another array with the relative displacements
 * for each chunk.
 */
void getDataDistribution (int n, int howMany, int **outFilasEnCPU, int **outDesplParaCPU) {
  int *filasEnCPU, *desplParaCPU;
  int filasPorCPU, filasPendientes, despl, i;
  
  /* Generate distribution of rows */
  filasEnCPU   = (int *)smalloc (howMany * sizeof (int));
  desplParaCPU = (int *)smalloc (howMany * sizeof (int));
  filasPorCPU  = n / howMany;
  filasPendientes = n - (filasPorCPU * howMany);
  despl = 0;
  for (i = 0; i < howMany; i++) {
    if (filasPendientes > 0) {
      filasEnCPU[i] = filasPorCPU + 1;
      filasPendientes--;
    } else {
      filasEnCPU[i] = filasPorCPU;
    }
    desplParaCPU[i] = despl;
    despl += filasEnCPU[i];
  }
  
  /* Return values */
  *outFilasEnCPU   = filasEnCPU;
  *outDesplParaCPU = desplParaCPU;
}
