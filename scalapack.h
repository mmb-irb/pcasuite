#ifndef SCALAPACK_H_
#define SCALAPACK_H_

/* Preprocessing of BLACS/ScaLAPACK symbols */
#ifdef UNDERSCORES_FORTRAN_
#define blacs_pinfo    blacs_pinfo__
#define blacs_get      blacs_get__
#define blacs_gridinit blacs_gridinit__
#define blacs_gridinfo blacs_gridinfo__
#define blacs_gridexit blacs_gridexit__
#define blacs_exit     blacs_exit__
#define igebs2d        igebs2d_
#define igebr2d        igebr2d_
#define dgebs2d        dgebs2d_
#define dgebr2d        dgebr2d_
#define dtrbs2d        dtrbs2d_
#define dtrbr2d        dtrbr2d_

#define numroc   numroc_
#define descinit descinit_
#define pdelset  pdelset_
#define pdelget  pdelget_
#define pdsyev   pdsyev_
#define pdsyevx  pdsyevx_
#endif


/* BLACS declarations */
void blacs_pinfo (int *iam, int *nprocs);
void blacs_get (int *ctxt, int *what, int *val);
void blacs_gridinit (int *ctxt, char *scope, int *nprow, int *npcol);
void blacs_gridinfo (int *ctxt, int *nprow, int *npcol, int *nprows, int *npcols);
void blacs_gridexit (int *ctxt);
void blacs_exit (int *cont);
void igebs2d (int *ctxt, char *scope, char *top, int *m, int *n, int *A, int *lda);
void igebr2d (int *ctxt, char *scope, char *top, int *m, int *n, int *A, int *lda, int *rsrc, int *csrc);
void dgebs2d (int *ctxt, char *scope, char *top, int *m, int *n, double *A, int *lda);
void dgebr2d (int *ctxt, char *scope, char *top, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc);
void dtrbs2d (int *ctxt, char *scope, char *top, char *uplo, char *diag, int *m, int *n, double *A, int *lda);
void dtrbr2d (int *ctxt, char *scope, char *top, char *uplo, char *diag, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc);

/* ScaLAPACK declarations */
int numroc (int *n, int *nb, int *roc, int *isrcproc, int *nprocs);
void descinit (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ctxt, int *lld, int *info);
void pdelset (double *A, int *ia, int *ja, int *descA, double *alpha);
void pdelget (char *scope, char *top, double *alpha, double *A, int *ia, int *ja, int *descA);
void pdsyev (char *jobz, char *uplo, int *n, double *A, int *ia, int *ja, int *descA, double *W, double *Z, int *iz, int *jz, int *descZ, double *work, int *lwork, int *info);
void pdsyevx (char *jobz, char *range, char *uplo, int *N, double *A, int *ia, int *ja, int *descA, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, int *nz, double *W, double *orfac, double *Z, int *iz, int *jz, int *descZ, double *work, int *lwork, int *iwork, int *liwork, int *ifail, int *iclustr, double *gap, int *info);

#endif /*SCALAPACK_H_*/
