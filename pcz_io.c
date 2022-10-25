#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "pcz_io.h"
#include "pczcompute.h"
#include "constants.h"
#include "utils.h"

/* This file contains all the routines that operates with the PCZ files. All
 * values obtained directly from the file are returned by routines in this
 * module. All values derived from the data contained in the file are computed
 * by routines in other modules, like pczcompute.h */

/* This method opens a PCZ-format file and reads the headers, filling the
 * PCZFILE structure with proper values */
PCZFILE *pczopen (char *fname) {
  PCZFILE *pczf;
  int nReaded, i;
  float *evals;

  /* Create and initialize the PCZFILE structure */
  pczf = (PCZFILE *)smalloc (sizeof (PCZFILE));
  strcpy (pczf->filename, fname);
  strcpy (pczf->title, "");
  strncpy (pczf->version, "????", 4);
  pczf->pczVersion     = 0;
  pczf->status         = 0;
  pczf->natoms         = 0;
  pczf->nvecs          = 0;
  pczf->nframes        = 0;
  pczf->vtot           = 0.0;
  pczf->vexplained     = 0.0;
  pczf->quality        = 0.0;
  pczf->dimensionality = 0;

  /* Open the file */
  pczf->file = fopen (pczf->filename, "rb");
  if (pczf->file == NULL) {
    perror ("Error: Cannot open input file");
    pczf->status = errno;
    return pczf;
  }

  /* Read key information */
  nReaded = 0;
  nReaded += fread ( pczf->version,        4,              1, pczf->file);
  if (strncmp (pczf->version, "PCZ", 3) != 0) {
    fprintf (stderr, "Error: Input file is not in PCZ format\n");
    pczf->status = -3;
    return pczf;
  }
  pczf->pczVersion = pczf->version[3] - '0';

  if (pczf->pczVersion == 4) {
	  nReaded += fread ( pczf->title,          80,             1, pczf->file);
	  nReaded += fread (&pczf->natoms,         sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->nframes,        sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->nvecs,          sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->vtot,           sizeof (float), 1, pczf->file);
	  nReaded += fread (&pczf->vexplained,     sizeof (float), 1, pczf->file);
	  nReaded += fread (&pczf->dimensionality, sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->rms,            sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->haveAtomNames,  sizeof (int),   1, pczf->file);
	  pczf->title[80] = (char)0;

	  if (am_big_endian ()) {
	    pczf->natoms         = swap_endianness (pczf->natoms);
	    pczf->nframes        = swap_endianness (pczf->nframes);
	    pczf->nvecs          = swap_endianness (pczf->nvecs);
	    pczf->vtot           = swap_endiannessf (pczf->vtot);
	    pczf->vexplained     = swap_endiannessf (pczf->vexplained);
	    pczf->dimensionality = swap_endianness (pczf->dimensionality);
	    pczf->rms            = swap_endianness (pczf->rms);
	    pczf->haveAtomNames  = swap_endianness (pczf->haveAtomNames);
	  }

	  if (pczf->haveAtomNames) {
	    pczf->atomNames = (PDBATOM *)smalloc (pczf->natoms * sizeof (PDBATOM));
	    for (i = 0; i < pczf->natoms; i++) {
	      fread (&pczf->atomNames[i].atomNumber,    sizeof (int),  1, pczf->file);
	      fread (&pczf->atomNames[i].atomName,      sizeof (char), 4, pczf->file);
	      fread (&pczf->atomNames[i].residueNumber, sizeof (int),  1, pczf->file);
	      fread (&pczf->atomNames[i].residueName,   sizeof (char), 3, pczf->file);
	      fread (&pczf->atomNames[i].chain,         sizeof (char), 1, pczf->file);
	      pczf->atomNames[i].atomName[4]    = (char)0;
	      pczf->atomNames[i].residueName[3] = (char)0;
	      rotateLeft(pczf->atomNames[i].atomName, 4);
	      trim (pczf->atomNames[i].atomName);
	      trim (pczf->atomNames[i].residueName);
	      if (am_big_endian ()) {
	        pczf->atomNames[i].atomNumber = swap_endianness (pczf->atomNames[i].atomNumber);
	        pczf->atomNames[i].residueNumber = swap_endianness (pczf->atomNames[i].residueNumber);
	      }
	    }
	  }

	  if (nReaded != 10) {
	    fprintf (stderr, "Error: Cannot read header data\n");
	    pczf->status = -1;
	    return pczf;
	  }
	  if (pczf->vtot == 0.0) {
	    fprintf (stderr, "Error: Total variance cannot be zero\n");
	    pczf->status = -2;
	    return pczf;
	  }
	  pczf->quality = (float)(100.0 * (pczf->vexplained / pczf->vtot));
	} else if (pczf->pczVersion == 3 || pczf->pczVersion == 2) {
	  nReaded += fread ( pczf->title,          80,             1, pczf->file);
	  nReaded += fread (&pczf->natoms,         sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->nframes,        sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->nvecs,          sizeof (int),   1, pczf->file);
	  nReaded += fread (&pczf->vtot,           sizeof (float), 1, pczf->file);
	  pczf->title[80] = (char)0;

	  if (am_big_endian ()) {
	    pczf->natoms         = swap_endianness (pczf->natoms);
	    pczf->nframes        = swap_endianness (pczf->nframes);
	    pczf->nvecs          = swap_endianness (pczf->nvecs);
	    pczf->vtot           = swap_endiannessf (pczf->vtot);
	  }

	  if (nReaded != 6) {
	    fprintf (stderr, "Error: Cannot read header data\n");
	    pczf->status = -1;
	    return pczf;
	  }
	  if (pczf->vtot == 0.0) {
	    fprintf (stderr, "Warning: Total variance cannot be zero\n");
	  }
	  if (pczf->pczVersion == 2) {
	    pczf->natoms /= 3;
	  }

	  evals = pczevals(pczf);
	  pczf->vexplained = 0.0;
	  for (i = 0; i < pczf->nvecs; i++)
	  	pczf->vexplained += evals[i];
    if (pczf->vtot != 0.0)
	    pczf->quality = (float)(100.0 * (pczf->vexplained / pczf->vtot));
    else
      pczf->quality = 0.0;
	} else {
	    fprintf (stderr, "Error: PCZ unexpected version format\n");
	    pczf->status = -6;
	    return pczf;
  }

  /*
  printf ("Title: %s\nnAtoms: %i\nnFrames: %i\nnVecs: %i\nQuality: %f%%\n",
    pczf->title, pczf->natoms, pczf->nframes, pczf->nvecs,
    pczf->quality);
  */

  return pczf;
}

/* This method closes the file and frees some data structures allocated
 * by the corresponding open function */
void pczclose (PCZFILE *pczf) {
  if (pczf->haveAtomNames)
    sfree (pczf->atomNames);
  sfree (pczf);
  fclose (pczf->file);
}

/* This method return the Eigen Values contained in the PCZ file */
float *pczevals (PCZFILE *pczf) {
  float *evals, *ptr;
  int i;

  evals = (float *)smalloc (pczf->nvecs * 4);
  ptr = evals;
  for (i = 0; i < pczf->nvecs; i++) {
    fseek (pczf->file, HDR_SIZE(pczf)+MEAN_SIZE(pczf)+PADDING(pczf)+i*(EIGEN_DATA_SIZE(pczf))+EVEC_SIZE(pczf), SEEK_SET);
    fread (ptr, 4, 1, pczf->file);
    ptr++;
  }

  if (am_big_endian ()) {
    for (i = 0; i < pczf->nvecs; i++) {
      evals[i] = swap_endiannessf (evals[i]);
    }
  }

  return evals;
}

/* This method returns the requested Eigen Vector contained in the PCZ file.
 * First vector is numberez zero (0) */
float *pczevec (PCZFILE *pczf, int nEvec) {
  float *evec;
  int i;

  evec = (float *)smalloc (EVEC_SIZE(pczf));
  fseek (pczf->file, HDR_SIZE(pczf)+MEAN_SIZE(pczf)+PADDING(pczf)+nEvec*EIGEN_DATA_SIZE(pczf), SEEK_SET);
  fread (evec, 4, pczf->natoms * 3, pczf->file);

  if (am_big_endian ()) {
    for (i = 0; i < pczf->natoms * 3; i++) {
      evec[i] = swap_endiannessf (evec[i]);
    }
  }

  return evec;
}

/* This method returns the projections of all the frames corresponding to
 * the requested Eigen Vector contained in the PCZ file */
float *pczproj (PCZFILE *pczf, int ip) {
  float *proj;
  int i;

  proj = (float *)smalloc (PROJ_SIZE(pczf));
  fseek (pczf->file, HDR_SIZE(pczf)+MEAN_SIZE(pczf)+PADDING(pczf)+ip*EIGEN_DATA_SIZE(pczf)+EVEC_SIZE(pczf)+EVAL_SIZE(pczf), SEEK_SET);
  fread (proj, 4, pczf->nframes, pczf->file);

  if (am_big_endian ()) {
    for (i = 0; i < pczf->nframes; i++) {
      proj[i] = swap_endiannessf (proj[i]);
    }
  }

  return proj;
}

/* This method returns the mean structure contained in the PCZ file */
float *pczavg (PCZFILE *pczf) {
  float *avg;
  int i;

  avg = (float *)smalloc (MEAN_SIZE(pczf));
  fseek (pczf->file, HDR_SIZE(pczf), SEEK_SET);
  fread (avg, 4, pczf->natoms * 3, pczf->file);

   if (am_big_endian ()) {
    for (i = 0; i < pczf->natoms * 3; i++) {
      avg[i] = swap_endiannessf (avg[i]);
    }
  }

 return avg;
}

/* This method returns the projections for all of the Eigen Vectors for the
 * requested frame */
float *pczscores (PCZFILE *pczf, int iframe) {
  float *scores;
  int i;

  scores = (float *)smalloc (pczf->nvecs   * sizeof (float));
  for (i = 0; i < pczf->nvecs; i++) {
    fseek (pczf->file, HDR_SIZE(pczf) + MEAN_SIZE(pczf) + PADDING(pczf) +
      i*EIGEN_DATA_SIZE(pczf) + EVEC_SIZE(pczf) + EVAL_SIZE(pczf) +
			(iframe * 4), SEEK_SET);
    fread (&scores[i], 4, 1, pczf->file);
  }

  if (am_big_endian ()) {
    for (i = 0; i < pczf->nvecs; i++) {
      scores[i] = swap_endiannessf (scores[i]);
    }
  }

 return scores;
}

/* This method writes a PCZ format file given its main components.
 * The method receives
 * fName: the name of the output file
 * pczf: the structure with some basic parameters filled
 * mean: a float vector with the mean structures
 * evals: a float vector with the eigenvalues
 * evecs: a float* vector with the eigenvectors
 * proj: a float* vector with the projections
 * The evecs and proj vectors are in C format, that means that evecs[i] is the
 * address to the beginning of the vector with the desired data and do not need
 * to be related to evecs[i+1].
 */
void pczwrite (char *fName, PCZFILE *pczf, float *mean, float *evals,
  float **evecs, float **proj) {
  FILE *fout;
  int natom3, i, k;
    
  natom3 = pczf->natoms * 3;
  fout = fopen (fName, "wb");
    
  if (am_big_endian ()) {
    pczf->natoms         = swap_endianness (pczf->natoms);
    pczf->nframes        = swap_endianness (pczf->nframes);
    pczf->nvecs          = swap_endianness (pczf->nvecs);
    pczf->vtot           = swap_endiannessf (pczf->vtot);
    pczf->vexplained     = swap_endiannessf (pczf->vexplained);
    pczf->dimensionality = swap_endianness (pczf->dimensionality);
    pczf->rms            = swap_endianness (pczf->rms);
    pczf->haveAtomNames  = swap_endianness (pczf->haveAtomNames);
  }
  
  /* Write header */
  fwrite (VERSION,               4,            1,   fout);
  fwrite (pczf->title,          80,            1,   fout);
  fwrite (&pczf->natoms,         1, sizeof (int),   fout);
  fwrite (&pczf->nframes,        1, sizeof (int),   fout);
  fwrite (&pczf->nvecs,          1, sizeof (int),   fout);
  fwrite (&pczf->vtot,           1, sizeof (float), fout);
  fwrite (&pczf->vexplained,     1, sizeof (float), fout);
  fwrite (&pczf->dimensionality, 1, sizeof (int),   fout);
  fwrite (&pczf->rms,            1, sizeof (int),   fout);
  fwrite (&pczf->haveAtomNames,  1, sizeof (int),   fout);
  if (am_big_endian ()) {
    pczf->natoms         = swap_endianness (pczf->natoms);
    pczf->nframes        = swap_endianness (pczf->nframes);
    pczf->nvecs          = swap_endianness (pczf->nvecs);
    pczf->vtot           = swap_endiannessf (pczf->vtot);
    pczf->vexplained     = swap_endiannessf (pczf->vexplained);
    pczf->dimensionality = swap_endianness (pczf->dimensionality);
    pczf->rms            = swap_endianness (pczf->rms);
    pczf->haveAtomNames  = swap_endianness (pczf->haveAtomNames);
  }

  if (pczf->haveAtomNames) {
    /* We write the atom IDs */
    for (i = 0; i < pczf->natoms; i++) {
      if (am_big_endian ()) {
        pczf->atomNames[i].atomNumber    = swap_endianness (pczf->atomNames[i].atomNumber);
        pczf->atomNames[i].residueNumber = swap_endianness (pczf->atomNames[i].residueNumber);
      }
      fwrite (&pczf->atomNames[i].atomNumber,     1, sizeof (int),  fout);
      fwrite ( pczf->atomNames[i].atomName,       4, sizeof (char), fout);
      fwrite (&pczf->atomNames[i].residueNumber,  1, sizeof (int),  fout);
      fwrite ( pczf->atomNames[i].residueName,    3, sizeof (char), fout);
      fwrite (&pczf->atomNames[i].chain,          1, sizeof (char), fout);
      if (am_big_endian ()) {
        pczf->atomNames[i].atomNumber    = swap_endianness (pczf->atomNames[i].atomNumber);
        pczf->atomNames[i].residueNumber = swap_endianness (pczf->atomNames[i].residueNumber);
      }
    }
  }

  /* Write the average structure */
  if (am_big_endian ()) {
    for (i = 0; i < natom3; i++) {
      mean[i] = swap_endiannessf (mean[i]);
    }
  }
  fwrite (mean, natom3, sizeof (float), fout);
  if (am_big_endian ()) {
    for (i = 0; i < natom3; i++) {
      mean[i] = swap_endiannessf (mean[i]);
    }
  }

  /* Write the eigenvectors, eigenvalues and projections */
  for (k = 0; k < pczf->nvecs; k++) {
    if (am_big_endian ()) {
      for (i = 0; i < natom3; i++) {
        evecs[k][i] = swap_endiannessf (evecs[k][i]);
      }
      evals[k] = swap_endiannessf (evals[k]);
      for (i = 0; i < pczf->nframes; i++) {
        proj[k][i] = swap_endiannessf (proj[k][i]);
      }
    }
    fwrite (evecs[k], natom3, sizeof (float), fout);
    fwrite (&evals[k], 1, sizeof (float), fout);
    fwrite (proj[k], pczf->nframes, sizeof (float), fout);
    if (am_big_endian ()) {
      for (i = 0; i < natom3; i++) {
        evecs[k][i] = swap_endiannessf (evecs[k][i]);
      }
      evals[k] = swap_endiannessf (evals[k]);
      for (i = 0; i < pczf->nframes; i++) {
        proj[k][i] = swap_endiannessf (proj[k][i]);
      }
    }
  }
  fclose (fout);
}
