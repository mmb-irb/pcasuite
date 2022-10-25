/*
 * PCAzip port from f90 source
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <sys/mman.h>
#include "pcazip.h"
#include "traj_io.h"
#include "pcz_io.h"
#include "trajutils.h"
#include "constants.h"
#include "utils.h"
#include "gaussianrms.h"
#include "scanner/aux_parser.h"
#ifdef MPI_VERSION_
#include "mpi_routines.h"
#include <mpi.h>
#endif


/* Global variable: controls verbosity */
int verbose;

int main (int argc, char **argv) {
  long int tiRead, tiFit, tiCovar, tiDiagonal, tiProjection, tiWriting;
  long int tfRead, tfFit, tfCovar, tfDiagonal, tfProjection, tfWriting;
  long int totalTime;
  FILE *fout;
  TRAJFILE xf;
  PCZFILE *pczf;
  PDBATOM *pdba, *npdba;
  CLParameters *clp;
  float  *x = NULL, *xt = NULL, vtot, vsum, *pn = NULL, ftmp;
  double *w = NULL, *z = NULL;
  double *cm;
  int natom3 = 0, nframes = 0, nvecs = 0, dimensionality;
  int i, j, k, ip1, haveAtomNames, nSelectedCandidates = 0;
  int *mask, nmask;
  char *title;

#ifdef MPI_VERSION_
  int whoAmI, howMany;
  double quality;
  char hostname[257];
#endif

  /* Do parameter recognition. All processes will have the command line
   * parameters available */
  clp = readCLParameters (argc, argv);

#ifdef MPI_VERSION_
  MPI_Init (&argc,&argv);

  MPI_Comm_size (MPI_COMM_WORLD, &howMany);
  MPI_Comm_rank (MPI_COMM_WORLD, &whoAmI);
  gethostname (hostname, 257);
  hostname[256] = (char)0;
  
  fprintf (stderr, "Process %i of %i (id: %i) running on %s\n", whoAmI + 1, howMany, whoAmI, hostname);
  MPI_Barrier (MPI_COMM_WORLD);
  
  if (whoAmI == 0) {
#endif

  /* Banner */
  if (clp->verbose)
    printf ("pcazip V%s\n", CODE_VERSION);    

  /* We initialize some variables and do parameter recognition */
  pdba          = NULL;
  haveAtomNames = FALSE;
  verbose = clp->verbose;
  
  /* Now we can check that we have correct parameters */
  checkCLP (clp);

  /* Reading PDB */
  if (clp->pdbfile != NULL) {
    pdba = (PDBATOM *)smalloc (clp->natom * sizeof (PDBATOM));
    pdbread (clp->pdbfile, clp->natom, pdba);
    haveAtomNames = TRUE;
  }

  /* Open input file and read header and basic information data */
  tiRead = getTime();
  natom3 = clp->natom * 3;

  trajopen (clp->infile, clp->natom, &xf);
  nframes = xf.nframes;
  title = xf.title;

  /* Mask out the unneeded atoms */
  mask = (int *)smalloc (natom3 * sizeof (int));
  if (clp->maskfile != NULL) {
    maskread (clp->maskfile, clp->natom, &nmask, mask);
    if (verbose) {
      fprintf (stderr, "Masking down to %i atoms\n", nmask);
    }
    clp->natom = nmask;
    natom3 = clp->natom * 3;
  } else {
    if (clp->maskString != NULL && haveAtomNames) {
      /* Mask processing */
      selectAtomMask (pdba, clp->natom, clp->maskString, &nmask, mask);
      if (clp->verbose) {
        fprintf (stderr, "Masking down to %i atoms\n", nmask);
      }
      clp->natom = nmask;
      natom3 = clp->natom * 3;
    } else {
      for (i = 0; i < natom3; i++) {
        mask[i] = TRUE;
      }
    }
  }

  /* Check that after the atom selection, the number of eigenvectors
   * requested are not too much */
  if (clp->nvecs > natom3 - 6) {
    fprintf (stderr, "Error: You are asking for %i eigenvectors, but you"
      " only have %i atoms. This is non-sense.\n", clp->nvecs, clp->natom);
    exit (1);
  }

  /* Use the atom selection to select the appropriate names */
  if (haveAtomNames) {
    npdba = (PDBATOM *)smalloc (clp->natom * sizeof (PDBATOM));
    packPdbAtom (npdba, pdba, mask, clp->natom);
    sfree (pdba);
    pdba = npdba;
  }

  /* Prepare the memory space for the matrixes. We ask for memory here and
   * not when needed because if we don't have all the memory needed, we can
   * avoid a computationally intensive calculus */
#ifndef MPI_VERSION_
  cm   = (double *)smalloc (natom3 * natom3 * sizeof (double));
#endif
  xt   = (float *)smalloc (xf.natoms * 3 * sizeof (float));

#ifdef USE_MMAP_
/**************************************/
/* MMap implementation to save memory */
/**************************************/
  int trajSize, error;
  FILE *tmpFile;
  
  /* File creation */
  trajSize = natom3 * nframes * sizeof (float);
  tmpFile = tmpfile(); /* This file will be automatically deleted */
  fseek (tmpFile, trajSize, SEEK_SET);
  error = fwrite (" ", 1, 1, tmpFile);
  fflush (tmpFile);
  if (error != 1) {
    fprintf (stderr, "Cannot allocate enough temporary disk space (%i bytes)\n", trajSize);
    exit (-3);
  }
  
  /* File mapping */
  x = (float *)mmap (NULL, trajSize, PROT_READ | PROT_WRITE, MAP_PRIVATE, fileno (tmpFile), 0);
#else
  x    = (float *)smalloc (natom3 * nframes * sizeof (float));
#endif

  /* Read all the snapshots, selecting only the desired atoms, and store on
   * memory for later fitting */
  for (i = 1; i <= nframes; i++) {
    trajsnap (&xf, i, xt);
    ip1 = natom3 * (i - 1);
    /*ip2 = ip1 + natom3;*/
    pack (x + ip1, xt, mask, xf.natoms * 3);
  }
  
  trajclose (&xf);
  if (clp->verbose)
    printf ("%i frames read\n", nframes);
  sfree (xt);
  tfRead = getTime();

  /* Fitting process: 1. fit to 1st frame, 2: calc avg, 3: fit to this */
  tiFit = getTime();
  if (clp->verbose)
    printf ("Fitting trajectory\n");
  xt = (float *)smalloc (natom3 * sizeof (*xt));
#ifdef MPI_VERSION_
  }
  if (clp->rms == STD_RMS) {
    /* Standard RMS */
    trajfit_mpi (x + natom3, x, clp->natom, nframes - 1);
    if (whoAmI == 0)
      trajavg (x, xt, clp->natom, nframes);
    trajfit_mpi (x, xt, clp->natom, nframes);
  } else if (clp->rms == GAUSSIAN_RMS) {
    /* Gaussian RMS */
    gaussianfit_mpi (x, x + natom3, clp->natom, nframes - 1);
    if (whoAmI == 0)
      trajavg (x, xt, clp->natom, nframes);
    gaussianfit_mpi (xt, x, clp->natom, nframes);
  }
  if (whoAmI == 0) {
#else
  if (clp->rms == STD_RMS) {
    /* Standard RMS */
    trajfit (x + natom3, x, clp->natom, nframes - 1);
    trajavg (x, xt, clp->natom, nframes);
    trajfit (x, xt, clp->natom, nframes);
  } else if (clp->rms == GAUSSIAN_RMS) {
    /* Gaussian RMS */
    gaussianfit (x, x + natom3, clp->natom, nframes - 1);
    trajavg (x, xt, clp->natom, nframes);
    gaussianfit (xt, x, clp->natom, nframes);
  }
#endif
  tfFit = getTime();

  /* Calculate covariance matrix */
  if (clp->verbose)
    printf ("Calculating covariance matrix\n");

  tiCovar = getTime();
#ifdef MPI_VERSION_
  }
#endif

  tiCovar = getTime();
#ifdef MPI_VERSION_
  meanCovar (x, nframes, natom3, xt, &cm);
#else
  calccovar (x, xt, cm, natom3, nframes);
#endif
  tfCovar = getTime();
  
#ifdef MPI_VERSION_
  if (whoAmI == 0) {
    quality = (double)clp->iqual;
#endif
  /* Start with algebra */

  /* If nvecs==0 we need to calc all eigenvectors to work out how
   * many are needed for each specified quality limit */
  if (clp->verbose && clp->nvecs == 0)
    printf ("Calculating eigenvalues and eigenvectors\n");
  else if (clp->verbose && clp->nvecs != 0)
    printf ("Calculating dimensionality\n");
#ifdef MPI_VERSION_
  }
#endif
  
  /* Eigenvalues calculus */
  tiDiagonal = getTime();
#ifdef MPI_VERSION_
    computeEigen (cm, natom3, &quality, &w, &z, &clp->nvecs, &vtot, &vsum, &dimensionality);

    if (whoAmI==0) {
      if (clp->verbose) {
        printf ("%8i eigenvectors has been computed\n", clp->nvecs);
        printf ("%8.4f%% of variance captured\n", (vsum * 100 / vtot));
      }
    }
#else
  w = computeEigenvalues (cm, natom3);
  
  /* Now calc the total variance */
  vtot = getTotalVariance (w, natom3);

  /* Now calc how many eigenvectors will be needed and the variance explained */
  if (clp->iqual != 0) {
    /* Selection based on the quality */
    nSelectedCandidates = selectEigenvaluesByQuality (w, natom3, (float)clp->iqual);
    if (clp->nvecs == 0) {
      /* Selection only based on the quality */
      clp->nvecs = nSelectedCandidates;
    } else {
      /* Selection based on quality AND eigenvectors requested.
       * We must have a minimum quality and a minimum number of eigenvectors */
      if (nSelectedCandidates > clp->nvecs) {
        clp->nvecs = nSelectedCandidates;
      }
    }
  }
  vsum = getVarianceForNEigenvalues (w, natom3, clp->nvecs);

  if (clp->verbose) {
    printf ("%8i eigenvectors will be used\n", clp->nvecs);
    printf ("%8.4f%% of variance captured\n", (vsum * 100 / vtot));
  }

  /* Dimensionality calculus */
  dimensionality = 0;
  while (w[natom3 - dimensionality - 1] >= 1.0) {
    dimensionality++;
  }

  /* Now the real thing: compute the eigenvectors needed for the explained
   * variance */
  if (clp->verbose) {
    printf ("Calculating eigenvectors\n");
  }

  z = computeEigenvectors (cm, natom3, clp->nvecs);
#endif
  tfDiagonal = getTime();

#ifdef MPI_VERSION_
  if (whoAmI == 0) {
#endif
  /* Free some expensive resources no longer needed */
  sfree (cm);

  /* Now we can do the projections - re-use work() for this */
  if (clp->verbose)
    printf ("Calculating projections\n");
  tiProjection = getTime();
  pn =  (float *)smalloc (clp->nvecs * nframes * sizeof (float));
  nvecs = clp->nvecs;
#ifdef MPI_VERSION_
  }
#endif
#ifdef MPI_VERSION_
  calcproj_mpi (natom3, nframes, nvecs, x, xt, z, pn);
#else
  calcproj (natom3, nframes, clp->nvecs, x, xt, z, pn);
#endif

#ifdef MPI_VERSION_
  if (whoAmI == 0) {
#endif
  tfProjection = getTime();

  /* Free some expensive resources no longer needed */
#ifdef USE_MMAP_
  munmap (x, trajSize);
  fclose (tmpFile);
#else
  sfree (x);
#endif
  
  /* Now write the compressed output file */
  tiWriting = getTime();
  if (clp->verbose)
    printf ("Writing compressed file\n");
  if (clp->formatted) {
    /* ASCII formatted output file */
    printf ("Escrivint fitxer formatejat %s\n", clp->outfile);
    fout = fopen (clp->outfile, "w");
    fprintf (fout, "%80s\n", title);
    fprintf (fout, "%8i%8i%8i%8.3f\n", clp->natom, nframes, clp->nvecs, vtot);
    for (i = 0; i < natom3; i++) {
      fprintf (fout, "%8.3f", xt[i]);
    }
    for (k = clp->nvecs - 1; k >= 0; k--) {
      fprintf (fout, "\n");
      for (i = k*natom3; i < (k+1)*natom3; i++) {
        ftmp = (float)z[i];
        fprintf (fout, "%8.3f", ftmp);
      }
      ftmp = (float)w[k];
      fprintf (fout, "\n%8.3f\n", ftmp);
      for (j = k*nframes; j < (k+1)*nframes; j++) {
        fprintf (fout, "%8.3f", pn[j]);
      }
      fprintf (fout, "\n");
    }
  } else {
    /* Binary output file */
    /* First we copy ans transform the matrices to reduce precision from
     * double to float. This transformation also lays the data in a more
     * logical way, having the most important eigenvector in the first row
     * of the matrix */
    float **adpEvecs = (float **)smalloc (clp->nvecs * sizeof (*adpEvecs));
    float **adpProj = (float **)smalloc (nframes * sizeof (*adpProj));
    float *adpEvals = (float *)smalloc (clp->nvecs * sizeof (*adpEvals));
    for (i = 0; i < clp->nvecs; i++) {
      adpEvecs[i] = (float *)smalloc (natom3 * sizeof (**adpEvecs));
      adpProj[i] = (float *)smalloc (nframes * sizeof (**adpProj));
      for (j = 0; j < natom3; j++) {
        adpEvecs[i][j] = (float)z[(clp->nvecs - 1 - i)*natom3+j];
      }
      adpEvals[i] = (float)w[natom3 - 1 - i];
      for (j = 0; j < nframes; j++) {
        adpProj[i][j] = (float)pn[(clp->nvecs - 1 - i)*nframes+j];
      }
    }
    pczf = (PCZFILE *)smalloc (sizeof (*pczf));
    strncpy (pczf->title, title, 81);
    strncpy (pczf->version, VERSION, 4);
    pczf->natoms = clp->natom;
    pczf->nvecs = clp->nvecs;
    pczf->nframes = nframes;
    pczf->vtot = vtot;
    pczf->vexplained = vsum;
    pczf->dimensionality = dimensionality;
    pczf->rms = clp->rms;
    pczf->haveAtomNames = haveAtomNames;
    pczf->atomNames = pdba;
    
    pczwrite (clp->outfile, pczf, xt, adpEvals, adpEvecs, adpProj);
    
    for (i = 0; i < clp->nvecs; i++) {
      sfree (adpEvecs[i]);
      sfree (adpProj[i]);
    }
    sfree (adpEvecs);
    sfree (adpProj);
    sfree (adpEvals);
    sfree (pczf);
  }
  tfWriting = getTime();
  
  if (clp->verbose) {
    totalTime = tfRead - tiRead + tfFit - tiFit + tfCovar - tiCovar +
      tfDiagonal - tiDiagonal + tfProjection - tiProjection + tfWriting - tiWriting;
    printf ("Timing summary (in microseconds):\n");
    printf ("Read data                 : %10li (%6.2f%%)\n", tfRead - tiRead, (float)(tfRead - tiRead)/(float)totalTime*100.0);
    printf ("Trajectory fitting        : %10li (%6.2f%%)\n", tfFit - tiFit, (float)(tfFit - tiFit)/(float)totalTime*100.0);
    printf ("Covariance matrix calculus: %10li (%6.2f%%)\n", tfCovar - tiCovar, (float)(tfCovar - tiCovar)/(float)totalTime*100.0);
    printf ("Eigenvector calculus      : %10li (%6.2f%%)\n", tfDiagonal - tiDiagonal, (float)(tfDiagonal - tiDiagonal)/(float)totalTime*100.0);
    printf ("Trajectory projection     : %10li (%6.2f%%)\n", tfProjection - tiProjection, (float)(tfProjection - tiProjection)/(float)totalTime*100.0);
    printf ("Write output              : %10li (%6.2f%%)\n", tfWriting - tiWriting, (float)(tfWriting - tiWriting)/(float)totalTime*100.0);
  }


  /* Resource freeing */
  sfree (pn);
  sfree (z);
  sfree (w);
  sfree (mask);
  sfree (xt);
  if (pdba != NULL)
    sfree (pdba);
  trajfree (&xf);
  freeCLParameters (clp);
#ifdef MPI_VERSION_  
  }
 MPI_Finalize ();
#endif

  return 0;
}

/*****************************************************************************/
/*                           End of main program                             */
/*****************************************************************************/

/* This method packs the mask into an array */
void pack (float *dest, float *source, int *mask, int length) {
  int i;

  for (i = 0; i < length; i++) {
    if (mask[i] != FALSE) {
      *dest = *source;
      dest++;
    }
    source++;
  }
}

/* This method packs the mask into an array */
void packPdbAtom (PDBATOM *dest, PDBATOM *source, int *mask, int length) {
  int i, copied;

  for (i = 0, copied = 0; copied < length; i++) {
    if (mask[i*3] != FALSE) {
      memcpy (dest, source, sizeof (PDBATOM));
      dest++;
      copied++;
    }
    source++;
  }
}

/* Method to read a mask from a PDB file */
void maskread (char *maskfile, int natom, int *nmask, int *mask) {
  FILE *mf;
  char line[82];
  char strAtomNumber[6];
  int atomNumber;

  *nmask = 0;
  memset (mask, 0, natom * 3 * sizeof (float));
  strAtomNumber[5] = (char)0;

  mf = fopen (maskfile, "r");
  fgets (line, 82, mf);
  while (!feof (mf)) {
    if (strncmp (line, "ATOM  ", 6) == 0) {
      (*nmask)++;
      strncpy (strAtomNumber, line+6, 5);
      atomNumber = atoi (strAtomNumber);
      mask[3*(atomNumber-1)]   = TRUE;
      mask[3*(atomNumber-1)+1] = TRUE;
      mask[3*(atomNumber-1)+2] = TRUE;
    }
    fgets (line, 82, mf);
  }
  fclose (mf);
}

/* Method to read a mask from a PDB file */
int atomsInPDB (char *pdbfile) {
  FILE *pdbf;
  char line[82];
  int nAtoms;

  nAtoms = 0;
  pdbf = fopen (pdbfile, "r");
  fgets (line, 82, pdbf);
  while (!feof (pdbf)) {
    if (strncmp (line, "ATOM  ", 6) == 0)
      nAtoms++;
    fgets (line, 82, pdbf);
  }
  fclose (pdbf);

  return nAtoms;
}

CLParameters *readCLParameters (int argc, char **argv) {
  static struct option long_options[] = {
    {"infile",       1, 0, 'i'},
    {"outfile",      1, 0, 'o'},
    {"natoms",       1, 0, 'n'},
    {"verbose",      0, 0, 'v'},
    {"mask",         1, 0, 'm'},
    {"maskString",   1, 0, 'M'},
    {"eigenvectors", 1, 0, 'e'},
    {"quality",      1, 0, 'q'},
    {"gaussian",     0, 0, 'g'},
    {"pdb",          1, 0, 'p'},
    {"formatted",    0, 0, 'f'},
    {"help",         0, 0, 'h'},
    {0,              0, 0, 0}
  };
  int c, option_index;
  char *strParams = "i:o:n:vm:M:e:q:gp:fh";
  CLParameters *clp;

  option_index = 0;
  clp = (CLParameters *)smalloc (sizeof (*clp));
  clp->nvecs         = 0;
  clp->iqual         = 0;
  clp->natom         = 0;
  clp->verbose       = FALSE;
  clp->formatted     = FALSE;
  clp->infile        = NULL;
  clp->outfile       = NULL;
  clp->maskfile      = NULL;
  clp->maskString    = NULL;
  clp->pdbfile       = NULL;
  clp->rms           = STD_RMS;

  /* Loop for recognizing all input parameters */
  for (c = getopt_long (argc, argv, strParams,
            long_options, &option_index);
       c != -1;
       c = getopt_long (argc, argv, strParams,
            long_options, &option_index)) {
    switch (c) {
    case 'i':
      clp->infile = copyString (optarg);
      break;
    case 'o':
      clp->outfile = copyString (optarg);
      break;
    case 'n':
      clp->natom = atoi (optarg);
      break;
    case 'v':
      clp->verbose = TRUE;
      break;
    case 'm':
      clp->maskfile = copyString (optarg);
      break;
    case 'M':
      clp->maskString = copyString (optarg);
      break;
    case 'e':
      clp->nvecs = atoi (optarg);
      break;
    case 'q':
      clp->iqual = atoi (optarg);
      break;
    case 'g':
      clp->rms = GAUSSIAN_RMS;
      break;
    case 'p':
      clp->pdbfile = copyString (optarg);
      break;
    case 'f':
      clp->formatted = TRUE;
      break;
    case 'h':
      helptext ();
      exit (0);
      break;
    default:
      printf ("Unrecognized option %c(%i)\n", c, c);
      exit (1);
      break;
    }
  }
  
  return clp;
}


void freeCLParameters (CLParameters *clp) {
  sfree (clp->infile);
  sfree (clp->outfile);
  sfree (clp->maskfile);
  sfree (clp->maskString);
  sfree (clp->pdbfile);
  sfree (clp);
}

/* Information method fot the user */
void helptext (void) {
  printf("Usage:\npcazip -i infile -o outfile -n natoms\n");
  printf("     [-v] [--mask maskfile] [-e nev] [-q qual] [--pdb pdbfile]\n");
  printf("Details:\n");
  printf("-i infile   | amber format trajectory input file\n");
  printf("-o outfile  | compressed output file\n");
  printf("-n natoms   | number of atoms in a snapshot\n");
  printf("-v          | verbose diagnostics\n");
  printf("-m mfile    | pdb format mask file - only atoms in this file\n");
  printf("            | will be included in the compression. Atom\n");
  printf("            | numbers (2nd field) must be correct!\n");
  printf("-M string   | ptraj-like format mask\n");
  printf("-e nev      | Include nev eigenvectors.\n");
  printf("-q qual     | Include enough eigenvectors to capture qual%%\n");
  printf("            | (1<=qual<=99) of the total variance\n");
  printf("-g          | Use a gaussian RMSd for fitting\n");
  printf("-p          | Use a PDB file for atom number and name extraction\n");
  printf("-h          | Shows this help\n");
  printf("Note: in the absence of any -q or -e option, default is to\n");
  printf("      include enough eiegnevectors to capture 90%% of variance\n");
}

void checkCLP (CLParameters *clp) {
  int pdbNAtoms;
  
  if (clp->infile == NULL || clp->outfile == NULL ||
    (clp->natom == 0 && clp->pdbfile == NULL)) {
    helptext ();
    fprintf (stderr, "Error:\nYou must specify md input file, "
       "compressed output file and number of atoms\n");
    exit (1);
  }

  if (!fexist (clp->infile)) {
    fprintf (stderr, "Error:\nCannot read trajectory file (%s)\n", clp->infile);
    exit (1);
  }
  if ((clp->natom == 0 && !fexist (clp->pdbfile))) {
    fprintf (stderr, "Error:\nCannot read PDB file (%s)\n", clp->pdbfile);
    exit (1);
  }
  if (clp->maskfile != NULL && !fexist (clp->maskfile)) {
    fprintf (stderr, "Error:\nCannot read mask file (%s)\n", clp->maskfile);
    exit (1);
  }

/*  
  if (clp->nvecs != 0 && clp->iqual != 0) {
    fprintf (stderr, "Error:\nYou can only ask for a number of vectors "
       "or a quality setting\n");
    exit (1);
  }
*/
  if (clp->nvecs == 0 && clp->iqual == 0) {
    clp->iqual = IDQ;
  }
/*
  if (clp->iqual != 0)
    clp->nvecs = 0;
*/
  if (clp->iqual < 0 || clp->iqual > 100) {
    fprintf (stderr, "Error:\n-q must be in range 1-99\n");
    exit (1);
  }

  /* Get number of atoms from PDB file */
  if (clp->pdbfile != NULL) {
    pdbNAtoms = atomsInPDB (clp->pdbfile);
    if (pdbNAtoms != clp->natom && clp->natom != 0) {
      fprintf (stderr, "Error: Atoms in PDB (%i) differ from number of atoms specified (%i)\n",
        pdbNAtoms, clp->natom);
      exit (1);
    }
    clp->natom = pdbNAtoms;
  }

  /* We show the user which parameters we readed */
  if (clp->verbose) {
    printf ("Parameters readed:\n");
    printf ("\tInput file: %s\n", clp->infile);
    printf ("\tOutput file: %s\n", clp->outfile);
    printf ("\tNumber of atoms: %i\n", clp->natom);
    if (clp->maskfile != NULL)
      printf ("\tMask file: %s\n", clp->maskfile);
    if (clp->maskString != NULL)
      printf ("\tMask: %s\n", clp->maskString);
    if (clp->pdbfile != NULL)
      printf ("\tPDB file: %s\n", clp->pdbfile);
    if (clp->nvecs != 0)
      printf ("\tNumber of vectors: %i\n", clp->nvecs);
    if (clp->iqual != 0)
      printf ("\tQuality setting: %i\n", clp->iqual);
    if (clp->rms != 0)
      printf ("\tGaussian RMS: %i\n", clp->rms);
  }
}

#ifndef MPI_VERSION_
/* Method for projection calculus
 * Parameters:
 * natom3  -> (Input)  Number of atoms *3
 * nframes -> (Input)  Number of frames in trajectory
 * nvecs   -> (Input)  Number of eigenvectors
 * x       -> (Input)  Trajectory
 * xt      -> (Input)  Mean structure
 * z       -> (Input)  Eigenvectors
 * pn      -> (Output) Projection
 */
void calcproj (int natom3, int nframes, int nvecs, float *x, float *xt,
  double *z, float *pn) {
  int i, j, k, l;

  for (k = 0; k < nframes; k++) {
    /* First substract average structure */
    for (i = 0; i < natom3; i++) {
      x[k * natom3 + i] -= xt[i];
    }

    /* Loop over eigenvectors */
    for (j = 0; j < nvecs; j++) {
      pn[j * nframes + k] = 0.0;
      for (l = 0; l < natom3; l++) {
        pn[j * nframes + k] += z[j * natom3 + l] * x[k * natom3 + l];
      }
    }
  }
}

/* Method for restoring the matrix diagonal */
void cmrestore (double *cm, double *cmd, int n) {
  int i;

  for (i = 0; i < n; i++) {
    cm[i*n+i] = cmd[i];
  }
}

/* Method for saving the matrix diagonal */
void cmsave (double *cm, double *cmd, int n) {
  int i;

  for (i = 0; i < n; i++) {
    cmd[i] = cm[i*n+i];
  }
}

/* This method computes the covariance matrix */
void calccovar (float *x, float *xt, double *cm, int natom3, int nframes) {
  double dx1, dx2, rframes;
  int i, j, k, kxnatom3, ixnatom3;
  
  #pragma omp parallel shared(cm, x, xt) private(k, i, j, kxnatom3, ixnatom3, dx1, dx2, rframes) firstprivate(nframes, natom3)
  {
  #pragma omp for schedule(static,100)
  for (i = 0; i < natom3 * natom3; i++) {
    cm[i] = 0.0;
  }
  
  for (k = 0; k < nframes; k++) {
    kxnatom3 = k * natom3;
    #pragma omp for schedule(static,10)
    for (i = 0; i < natom3; i++) {
      ixnatom3 = i * natom3;
      dx1 = x[kxnatom3+i] - xt[i];
      for (j = i; j < natom3; j++) {
        dx2 = x[kxnatom3+j] - xt[j];
        cm[ixnatom3+j] += dx1 * dx2;
      }
    }
  }

  /* Reduce to average and fill the other triangle */
  rframes = 1.0 / nframes;
  #pragma omp for schedule(static,10)
  for (i = 0; i < natom3; i++) {
    for (j = i; j < natom3; j++) {
      cm[i*natom3+j] *= rframes;
      cm[j*natom3+i] = cm[i*natom3+j];
    }
  }
  } /* OMP pragma end */
}

double *computeEigenvalues (double *cm, int natom3) {
  char V, I, U, N, A;
  double vl, vu, abstol, *work, *w, *z, *cmd;
  int m, il, iu, lwork, liwork, info, *iwork, *isuppz;
  
  V = 'V';
  I = 'I';
  U = 'U';
  N = 'N';
  A = 'A';
  vl = 0.0;
  vu = 0.0;
  il = 1;
  iu = natom3;
  abstol = 0.0;
  lwork = -1;
  liwork = -1;
  work   = (double *)smalloc (sizeof (double));
  iwork  = (int *)smalloc (sizeof (int));
  isuppz = (int *)smalloc (2 * natom3 * sizeof (int));
  w      = (double *)smalloc (natom3 * sizeof (double));
  z      = (double *)smalloc (sizeof (double));
  cmd    = (double *)smalloc (natom3 * sizeof (double));

  dsyevr_(&V, &I, &U, &natom3, cm, &natom3, &vl, &vu, &il, &iu, &abstol,
    &m, w, z, &natom3, isuppz, work, &lwork, iwork, &liwork, &info);
  if (info != 0) {
    fprintf (stderr, "DSYEVR returns with info=%i\n", info);
    exit (1);
  }
  
  lwork = (int)*work;
  liwork = *iwork;
  
  sfree (work);
  sfree (iwork);
  work = (double *)smalloc (lwork * sizeof (double));
  iwork = (int *)smalloc (liwork * sizeof (int));

  /* Copy the leading diagonal of cm into cmd, as it will be destroyed
   * in the first call to dsyevr_ */
  cmsave (cm, cmd, natom3);

  /* Compute the eigenvalues */
  vl = 0.0;
  vu = 0.0;
  il = 1;
  iu = natom3;
  dsyevr_(&N, &A, &U, &natom3, cm, &natom3, &vl, &vu, &il, &iu, &abstol,
    &m, w, z, &natom3, isuppz, work, &lwork, iwork, &liwork, &info);
  if (info != 0) {
    fprintf (stderr, "DSYEVR returns with info=%i\n", info);
    exit (1);
  }

  /* Restore the diagonal */
  cmrestore (cm, cmd, natom3);
  
  /* Resource freeing */
  sfree (work);
  sfree (iwork);
  sfree (isuppz);
  sfree (z);
  sfree (cmd);
  
  return w;
}

double *computeEigenvectors (double *cm, int natom3, int nvecs) {
  char V, I, N, A, L;
  double vl, vu, abstol, *work, *w, *z;
  int m, il, iu, lwork, liwork, info, *iwork, *isuppz;
  
  V = 'V';
  I = 'I';
  N = 'N';
  A = 'A';
  L = 'L';
  vl = 0.0;
  vu = 0.0;
  il = 1;
  iu = natom3;
  abstol = 0.0;
  lwork = -1;
  liwork = -1;
  work   = (double *)smalloc (sizeof (double));
  iwork  = (int *)smalloc (sizeof (int));
  isuppz = (int *)smalloc (2 * natom3 * sizeof (int));
  w      = (double *)smalloc (natom3 * sizeof (double));
  z      = (double *)smalloc (nvecs * natom3 * sizeof (double));
  
  dsyevr_(&V, &I, &L, &natom3, cm, &natom3, &vl, &vu, &il, &iu, &abstol,
    &m, w, z, &natom3, isuppz, work, &lwork, iwork, &liwork, &info);
  if (info != 0) {
    fprintf (stderr, "DSYEVR returns with info=%i\n", info);
    exit (1);
  }
  
  lwork = (int)*work;
  liwork = *iwork;
  
  sfree (work);
  sfree (iwork);
  work = (double *)smalloc (lwork * sizeof (double));
  iwork = (int *)smalloc (liwork * sizeof (int));

  /* Compute the eigenvectors */
  il = natom3 - nvecs + 1;
  iu = natom3;

  /* Note we use the lower triangle of cm this time... */
  dsyevr_(&V, &I, &L, &natom3, cm, &natom3, &vl, &vu, &il, &iu, &abstol,
    &m, w, z, &natom3, isuppz, work, &lwork, iwork, &liwork, &info);
  if (info != 0) {
    fprintf (stderr, "DSYEVR returns with info=%i\n", info);
    exit (1);
  }

  /* Resource freeing */
  sfree (work);
  sfree (iwork);
  sfree (isuppz);
  sfree (w);
  
  return z;
}
#endif

double getQualityForNEigenvalues (double *w, int length, int nEigenvalues) {
  double vtot, accVar, quality;
  
  vtot = getTotalVariance (w, length);
  accVar = getVarianceForNEigenvalues (w, length, nEigenvalues);
  quality = accVar * 100.0 / vtot;
  
  return quality;
}

double getVarianceForNEigenvalues (double *w, int length, int nEigenvalues) {
  double accVar;
  int i;
  
  accVar = 0.0;
  for (i = 0; i < nEigenvalues; i++) {
    accVar += w[length - 1 - i];
  }

  return accVar;
} 

int selectEigenvaluesByQuality (double *w, int length, float quality) {
  double vtot, accVar, epsilon;
  int i, selectedVecs;
  
  accVar = 0.0;
  vtot = getTotalVariance (w, length);
  epsilon = vtot * 0.0000000001;
  i = length - 1;
  while ((i >= 0) &&  /* Index check */ 
    (w[i] > epsilon) && /* Value check */
    ((accVar * 100.0 / vtot) < quality) /* Condition check */) {
    accVar += w[i];
    i--;
  }

  selectedVecs = length - 1 - i;
  
  return selectedVecs;
} 

double getTotalVariance (double *w, int length) {
  double vtot;
  int i;
  
  vtot = 0.0;
  for (i = 0; i < length; i++) {
    vtot += w[i];
  }
  
  return vtot;
}
