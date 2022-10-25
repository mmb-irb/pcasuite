/*
 * PCAunzip port from f90 source
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include "pcaunzip.h"
#include "utils.h"
#include "constants.h"
#include "pcz_io.h"
#include "pczcompute.h"
#include "scanner/parser.h"

int verbose;

/* Program to regenerate the original trajectory (with the generated compression
 * error) from the compressed file */
int main (int argc, char **argv) {
  int i, c, pdb, nAtoms, *atomList;
  char *infile, *outfile, *maskString;
  float *snapshot, *selection;
  FILE *fout;
  PCZFILE *pczf;
  int option_index = 0;
  static struct option long_options[] = {
    {"infile",  1, 0, 'i'},
    {"outfile", 1, 0, 'o'},
    {"pdb",     0, 0, 'p'},
    {"mask",    1, 0, 'M'},
    {"verbose", 0, 0, 'v'},
    {"help",    0, 0, 'h'},
    {0,         0, 0, 0}
  };

  /* We initialize some variables in order to do parameter recognition */
  infile     = NULL;
  outfile    = NULL;
  maskString = NULL;
  verbose    = FALSE;
  pdb        = FALSE;

  /* Loop for recognizing input parameters */
  for (c = getopt_long (argc, argv, "i:o:pM:vh",
      long_options, &option_index);
       c != -1;
       c = getopt_long (argc, argv, "i:o:pM:vh",
      long_options, &option_index)) {
    switch (c) {
    case 'i':
      infile = copyString (optarg);
      break;
    case 'o':
      outfile = copyString (optarg);
      break;
    case 'p':
      pdb = TRUE;
      break;
    case 'M':
      maskString = copyString (optarg);
      break;
    case 'v':
      verbose = TRUE;
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

  /* Banner */
  if (verbose)
    printf ("pcaunzip V%s\n", CODE_VERSION);  
  
  if (optind < argc) {
    printf ("There are unrecognized options in the command line: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
    printf ("\n");
  }

  /* Now we can check that we have correct parameters */
  if (infile == NULL) {
    helptext ();
    fprintf (stderr, "Error:\nYou must specify pcz input file\n");
    exit (1);
  } else {
    pczf = pczopen (infile);
    if (pczf->status != 0) {
      fprintf (stderr, "Error with input file, pczopen returns: %i\n",
         pczf->status);
      exit (1);
    }
  }

  if (outfile != NULL) {
    fout = fopen (outfile, "w");
  } else {
    fout = stdout;
  }

  /* We show the user which parameters we readed */
  if (verbose) {
    printf ("Parameters readed:\n");
    printf ("\tInput file: %s\n", infile);
    if (outfile != NULL)
      printf ("\tOutput file: %s\n", outfile);
    else
      printf ("\tOutput goes to standard output\n");
    if (maskString != NULL)
      printf ("\tMask: %s\n", maskString);
  }



  /* We process the mask in order to get an atom list */
  /* First we check that we have all the information necessary to do the
   * calculus */
  if (maskString != NULL) {
    if (!pczf->haveAtomNames) {
      fprintf (stderr, "Warning: No atom names are present, computing data with all atoms\n");
      atomList = NULL;
      nAtoms = pczf->natoms;
    } else {
      atomList = selectAtoms (pczf, maskString, &nAtoms);

      if (nAtoms == 0) {
        fprintf (stderr, "Error: There are no atoms in selection\n");
        return -2;
      }
      fprintf (stderr, "There are %i atoms in selection\n", nAtoms);
    }
  } else {
      nAtoms = pczf->natoms;
      atomList = NULL;
      fprintf (stderr, "Warning: No selection made, computing data with all atoms\n");
  }

  /* Finally, write the trajectory */
  if (pdb == FALSE) {
    /* Raw output format */
    fprintf (fout, "%-80.80s\n", trim(pczf->title));
    for (i = 0; i < pczf->nframes; i++) {
      snapshot = pczsnap (pczf, i);
      if (snapshot == NULL)
        return -1;
      if (atomList == NULL) {
        writeVector (snapshot, pczf->natoms * 3, fout);
      } else {
        selection = packAtoms (snapshot, atomList, nAtoms);
        writeVector (selection, nAtoms * 3, fout);
        sfree (selection);
      }
      sfree (snapshot);
    }
  } else {
    /* PDB output format */
    fprintf (fout, "TITLE     %-60.60s\n", trim(pczf->title));
    for (i = 0; i < pczf->nframes; i++) {
      fprintf (fout, "MODEL     %4i\n", i + 1);
      snapshot = pczsnap (pczf, i);
      if (snapshot == NULL)
        return -1;
      if (atomList == NULL) {
        dumppdb (pczf, fout, snapshot, NULL, NULL, pczf->natoms);
      } else {
        selection = packAtoms (snapshot, atomList, nAtoms);
        dumppdb (pczf, fout, selection, NULL, NULL, nAtoms);
        sfree (selection);
      }
      sfree (snapshot);
      fprintf (fout, "ENDMDL\n");
    }
  }

  pczclose (pczf);
  if (outfile != NULL)
    fclose (fout);

  /* Resource freeing */
  sfree (infile);
  if (outfile != NULL)
    sfree (outfile);

  return 0;
}

/* Information method fot the user */
void helptext (void) {
  printf("Usage:\npcaunzip -i infile [-o outfile] [--pdb] [--verbose] [--help]\n");
}
