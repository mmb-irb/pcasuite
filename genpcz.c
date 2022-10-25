#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "genpcz.h"
#include "constants.h"
#include "pcz_io.h"
#include "utils.h"

/*
 * The main use of this program is to generate a PCZ format file from
 * other text files containing the needed data in text form.
 * The data needed for a PCZ file to be generated is:
 *   -Atom and residue names and numbers
 *   -Mean structure
 *   -Eigenvalues
 *   -Eigenvectors
 *   -Projections
 * So the command line can be similar to:
 * $ ./genpcz -o res.pcz -m mean.pdb -e eigen.data -a evals.dat
 */
int main (int argc, char **argv) {
  float *evals, **evecs, *mean, **proj;
  int nEvals, nEvecs, nEvecsInProj;
  PCZFILE pczf;
  CLParameters *clp;
  int i;
  
  /* Parameter processing */
  clp = readCLParameters (argc, argv);
  
  /* Banner */
  if (clp->verbose)
    printf ("genpcz V%s\n", CODE_VERSION);  
  
  /* Eigenvalue load */
  evals = readEvals (clp->evalsFile, &nEvals);
  if (clp->verbose)
    printf ("%i atoms found\n", nEvals / 3);
  
  /* Eigenvector and mean load */
  evecs = readEvecs (clp->evecsFile, &nEvecs, &mean);
  if (clp->verbose) 
    printf ("%i eigenvectors found\n", nEvecs);

  /* Projections load */
  proj = readProj (clp->projFile, &nEvecsInProj, &pczf.nframes);
  if (clp->verbose) {
    printf ("%i vectors projected found\n", nEvecsInProj);
    printf ("%i frames found\n", pczf.nframes);
  }
  
  /* Variance calculus */
  pczf.nvecs = MIN (nEvecs, nEvecsInProj);
  pczf.vtot = 0.0;
  for (i = 0; i < nEvals; i++) {
    pczf.vtot += evals[i];
  }
  pczf.vexplained = 0.0;
  for (i = 0; i < pczf.nvecs; i++) {
    pczf.vexplained += evals[i];
  }
  pczf.quality = pczf.vexplained * 100 / pczf.vtot;
  if (clp->verbose) {
    printf ("Using %i eigenvectors\n", pczf.nvecs);
    printf ("Total variance: %f -- Explained variance: %f\n", pczf.vtot, pczf.vexplained);
    printf ("Quality: %f\n", pczf.quality);
  }
  
  /* Pending header data initialization */
  strncpy (pczf.title, clp->title, 81);
  pczf.title[80] = (char)0;
  strncpy (pczf.version, VERSION, 4);
  pczf.natoms = nEvals / 3;
  pczf.dimensionality = clp->dimensionality;
  pczf.rms = clp->rmsType;

  /* PDB atom name retrieval */
  pczf.atomNames = (PDBATOM *)smalloc (pczf.natoms * sizeof (PDBATOM));
  pdbread (clp->pdbFile, pczf.natoms, pczf.atomNames);
  pczf.haveAtomNames = TRUE;
  if (clp->verbose)
    printf ("Atom names readed\n");
  
  pczwrite (clp->pczFile, &pczf, mean, evals, evecs, proj);
  if (clp->verbose)
    printf ("Output file written\n");
  
  /* Resource freeing */
  sfree (evals);
  for (i = 0; i < nEvecs; i++)
    sfree (evecs[i]);
  sfree (evecs);
  sfree (mean);
  for (i = 0; i < nEvecsInProj; i++)
    sfree (proj[i]);
  sfree (proj);
  sfree (pczf.atomNames);
  freeCLParameters (clp);
  
  return 0;
}

float *readEvals (char *evalsFile, int *nEvals) {
  float *evals;
  FILE *fi;
  char *line;
  int i, idx;
  
  fi = fopen (evalsFile, "r");
  if (fi == NULL) {
    fprintf (stderr, "Error: Cannot open eigenvalues file (%s)\n", evalsFile);
    exit (1);
  }
  
  line = (char *)smalloc (MAX_LINE_LENGTH * sizeof (*line));
  line = fgets (line, MAX_LINE_LENGTH, fi);
  fscanf (fi, "%i %i\n", nEvals, nEvals);
  if (feof (fi)) {
    fprintf (stderr, "Error: Eigenvalues file empty (%s)\n", evalsFile);
    exit (1);
  }
  evals = (float *)smalloc (*nEvals * sizeof (*evals));
  for (i = 0; i < *nEvals; i++) {
    fscanf (fi, " ****\n");
    fscanf (fi, "%i %f", &idx, &evals[i]);
  }
  
  free (line);
  fclose (fi);
  
  return evals;
}

float **readEvecs (char *evecsFile, int *nEvecs, float **mean){
  float **evecs;
  FILE *fi;
  char *line;
  int nAtomsx3;
  int i, j, tmp;
  float ftmp;

  fi = fopen (evecsFile, "r");
  if (fi == NULL) {
    fprintf (stderr, "Error: Cannot open eigenvectors file (%s)\n", evecsFile);
    exit (1);
  }
 
  line = (char *)smalloc (MAX_LINE_LENGTH * sizeof (*line));
  line = fgets (line, MAX_LINE_LENGTH, fi);
  fscanf (fi, "%i %i Contains %i eigenvectors\n", &nAtomsx3, &nAtomsx3, nEvecs);
  if (feof (fi)) {
    fprintf (stderr, "Error: Eigenvectors file empty (%s)\n", evecsFile);
    exit (1);
  }
   
  /* Mean reading */
  *mean = (float *)smalloc (nAtomsx3 * sizeof (**mean));
  for (i = 0; i < nAtomsx3; i++) {
    fscanf (fi, "%f", &(*mean)[i]);
  }

  evecs = (float **)smalloc (*nEvecs * sizeof (*evecs));
  for (i = 0; i < *nEvecs; i++) {
    evecs[i] = (float *)smalloc (nAtomsx3 * sizeof (*evecs[i]));
    fscanf (fi, " ****\n"); 
    fscanf (fi, "%i %f\n", &tmp, &ftmp); 
    for (j = 0; j < nAtomsx3; j++) {
      fscanf (fi, "%f", &evecs[i][j]);
    }
  }
  
  sfree (line);
  fclose (fi);
  
  return evecs;
}

float **readProj (char *projFile, int *nEvecs, int *nFrames) {
  float **proj, ftmp;
  FILE *fi;
  int i, j, tmp;

  fi = fopen (projFile, "r");
  if (fi == NULL) {
    fprintf (stderr, "Error: Cannot open projections file (%s)\n", projFile);
    exit (1);
  }
  
  /* In the first place, we find how many modes we have projected */
  fscanf (fi, "Projection of snapshots onto modes\n");
  fscanf (fi, " Snapshot ");
  if (feof (fi)) {
    fprintf (stderr, "Error: Projections file empty (%s)\n", projFile);
    exit (1);
  }
  while (fscanf (fi, " Mode%i ", nEvecs) != 0)
    ;
  
  /* In second place we find how many frames has the file */
  *nFrames = 0;
  fscanf (fi, "%i", &tmp);
  while (!feof (fi)) {
    for (i = 0; i < *nEvecs; i++) {
      fscanf (fi, "%f", &ftmp);
    }
    (*nFrames)++;
    fscanf (fi, "%i", &tmp);
  }

  /* One we have the number of frames, we restart reading the file and
   * loading the values into memory */
  rewind (fi);
  fscanf (fi, "Projection of snapshots onto modes\n");
  fscanf (fi, " Snapshot ");
  while (fscanf (fi, " Mode%i ", nEvecs) != 0)
    ;
  
  /* Memory reservation for the projections */
  proj = (float **)smalloc (*nEvecs * sizeof (*proj));
  for (i = 0; i < *nEvecs; i++) {
    proj[i] = (float *)smalloc (*nFrames * sizeof (**proj));
  }
  
  /* Data reading */
  for (j = 0; j < *nFrames; j++) {
    fscanf (fi, "%i", &tmp);
    for (i = 0; i < *nEvecs; i++) {
      fscanf (fi, "%f", &proj[i][j]);
    }
  }
  
  fclose (fi);
  
  return proj;
}

CLParameters *readCLParameters (int argc, char **argv) {
  static struct option long_options[] = {
    {"evalsfile",       1, 0, 'a'},
    {"evecsfile",       1, 0, 'e'},
    {"projfile",        1, 0, 'j'},
    {"pdbfile",         1, 0, 'p'},
    {"pczfile",         1, 0, 'o'},
    {"title",           1, 0, 't'},
    {"dimensionality",  1, 0, 'd'},
    {"rms",             1, 0, 'r'},
    {"verbose",         0, 0, 'v'},
    {"help",            0, 0, 'h'},
    {0,                 0, 0, 0}
  };
  int c, option_index;
  char *strParams = "a:e:j:p:o:vh";
  CLParameters *clp;

  option_index = 0;
  clp = (CLParameters *)smalloc (sizeof (*clp));
  clp->evalsFile      = NULL;
  clp->evecsFile      = NULL;
  clp->projFile       = NULL;
  clp->pdbFile        = NULL;
  clp->pczFile        = NULL;
  clp->rmsType        = NONE_RMS;
  clp->dimensionality = -1;
  clp->title          = copyString ("genpcz generated file");
  clp->verbose        = FALSE;
  
  /* Loop for recognizing all input parameters */
  for (c = getopt_long (argc, argv, strParams,
              long_options, &option_index);
       c != -1;
       c = getopt_long (argc, argv, strParams,
              long_options, &option_index)) {
    switch (c) {
    case 'a':
      clp->evalsFile = copyString (optarg);
      break;
    case 'e':
      clp->evecsFile = copyString (optarg);
      break;
    case 'j':
      clp->projFile = copyString (optarg);
      break;
    case 'p':
      clp->pdbFile = copyString (optarg);
      break;
    case 'o':
      clp->pczFile = copyString (optarg);
      break;
    case 't':
      sfree (clp->title);
      clp->title = copyString (optarg);
      break;
    case 'd':
      clp->dimensionality = atoi (optarg);
      break;
    case 'r':
      if (strcmp(optarg, "std") == 0) {
        clp->rmsType = STD_RMS;
      } else if (strcmp(optarg, "gaussian") == 0) {
        clp->rmsType = GAUSSIAN_RMS;
      }
      break;
    case 'v':
      clp->verbose = TRUE;
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
  sfree (clp->evalsFile);
  sfree (clp->evecsFile);
  sfree (clp->projFile);
  sfree (clp->pdbFile);
  sfree (clp->pczFile);
  sfree (clp->title);
}

void helptext (void) {
  printf ("genpcz v0.01\n");
  printf ("  use: ./genpcz <parameters>\n");
  printf ("  Parameters:\n");
  printf ("  --evalsfile <name>     : Gives the file with the eigenvalues\n");
  printf ("  --evecsfile <name>     : Gives the file with the eigenvectors\n");
  printf ("  --projfile <name>      : Gives the file with the projections\n");
  printf ("  --pdbfile <name>       : Gives the file with the names of the atoms\n");
  printf ("  --pczfile <name>       : Gives the file that will contain the output\n");
  printf ("  --title <name>         : Gives the title for the contents\n");
  printf ("  --dimensionality <name>: Gives the dimensionality of the data\n");
  printf ("  --rms <name>           : Gives the RMSd type used [std|gaussian|none]\n");
  printf ("  --verbose              : Show progress information\n");
  printf ("  --help                 : Show this help message\n");
}
