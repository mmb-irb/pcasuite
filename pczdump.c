/*
 * PCZdump port from f90 source + new additions
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include "pczdump.h"
#include "utils.h"
#include "constants.h"
#include "pcz_io.h"
#include "pczcompute.h"
#include "mahalanobis.h"
#include "clusterlist.h"
#include "scanner/parser.h"

int verbose;
float threshold;

/* This is just the main program loop */
int main (int argc, char **argv) {
  int nAtoms, *atomList;
  FILE *fout;
  PCZFILE *pczf;
  CLParameters *clp;

  clp = readCLParameters (argc, argv);
  verbose = clp->verbose;
  
  /* Banner */
  if (clp->verbose)
    printf ("pczdump V%s\n", CODE_VERSION);  
  
  if (clp->infile == NULL) {
    helptext ();
    fprintf (stderr, "Error:\nYou must specify pcz input file\n");
    exit (1);
  } else {
    pczf = pczopen (clp->infile);
    if (pczf->status != 0) {
      fprintf (stderr, "Error with input file, pczopen returns: %i\n",
         pczf->status);
      exit (1);
    }
  }

  /* Command line parameter validation */
  checkCLP (clp, pczf);

  if (clp->outfile != NULL) {
    fout = fopen (clp->outfile, "w");
  } else {
    fout = stdout;
  }

  /* We process the mask in order to get an atom list */
  /* First we check that we have all the information necessary to do the
   * calculus */
  if (clp->maskString != NULL) {
    if (!pczf->haveAtomNames) {
      fprintf (stderr, "Warning: No atom names are present, computing data with all atoms\n");
      atomList = NULL;
      nAtoms = pczf->natoms;
    } else {
      atomList = selectAtoms (pczf, clp->maskString, &nAtoms);

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

  /* Finally, we call the requested routine */
  if (clp->averageRequested) {
    if (clp->pdbFormatRequested) {
      dumpavgpdb(pczf, fout, atomList, nAtoms/*, pdbline, nr*/);
    } else {
      dumpavgx (pczf, fout, atomList, nAtoms);
    }
  } else if (clp->evalsRequested != FALSE) {
    dumpevals (pczf, fout);
  } else if (clp->infoRequested != FALSE) {
    dumpinfo (pczf, fout);
  } else if (clp->nEvec != 0) {
    if (clp->pdbFormatRequested) {
      dumpevecpdb(pczf, fout, clp->nEvec - 1, atomList, nAtoms);
    } else {
      dumpevec (pczf, fout, clp->nEvec - 1, atomList, nAtoms);
    }
  } else if (clp->nProj != 0) {
    dumpproj (pczf, fout, clp->nProj - 1);
  } else if (clp->rmsRequested != FALSE) {
    dumprms (pczf  , fout, clp->nRms - 1);
  } else if (clp->lindemannRequested != FALSE) {
    dumplindemann (pczf, fout, atomList, nAtoms);
  } else if (clp->forceCteRequested != FALSE) {
    dumpforcecte (pczf, fout, clp->temperature);
  } else if (clp->collectivityRequested != FALSE) {
    dumpcollectivity (pczf, fout, clp->nEvecCollectivity - 1);
  } else if (clp->hingeRequested != FALSE) {
    dumphinge (pczf, fout, clp->nEvecHinge - 1);
  } else if (clp->stiffRequested != FALSE) {
    dumpstiffness (pczf, fout, clp->nEvecStiff - 1, nAtoms, atomList, clp->temperature);
  } else if (clp->mahalanobisRequested != FALSE) {
    dumpmahalanobis (pczf, fout, clp->nev, clp->xfile, atomList, nAtoms);
  } else if (clp->flucRequested != 0) {
    if (clp->pdbFormatRequested) {
      dumpflucpdb (pczf, fout, clp->nFluc - 1, clp->bfactorsRequested, atomList, nAtoms);
    } else {
      dumpfluc (pczf, fout, clp->nFluc - 1,clp-> bfactorsRequested, atomList, nAtoms);
    }
  } else if (clp->nAnim != 0) {
    if (clp->pdbFormatRequested) {
      dumpanimpdb (pczf, fout, clp->nAnim - 1, atomList, nAtoms/*, pdbline, nr*/);
    } else {
      dumpanimx (pczf, fout, clp->nAnim - 1, atomList, nAtoms);
    }
  }

  if (clp->outfile != NULL)
    fclose (fout);
  pczclose (pczf);

  /* Resource freeing */
  if (clp->maskString != NULL) {
    sfree (atomList);
  }
  freeCLParameters (clp);

  return 0;
}

/*****************************************************************************/
/*                           End of main program                             */
/*****************************************************************************/

/* This method dumps the eigenvectors to a file */
void dumpevec (PCZFILE *pczf, FILE *fout, int nEvec, int *atomList, int nAtoms) {
  float *evecs, *selection;

  evecs = pczevec (pczf, nEvec);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packAtoms (evecs, atomList, nAtoms);
    sfree (evecs);
    evecs = selection;
  }

  writeVector (evecs, nAtoms * 3, fout);
  sfree (evecs);
}

/* This method dumps the eigenvalues to a file */
void dumpevals (PCZFILE *pczf, FILE *fout) {
  float *evals;

  evals = pczevals (pczf);
  writeColumnVector (evals, pczf->nvecs, fout);
  sfree (evals);
}

/* This method dumps the projections to a file */
void dumpproj (PCZFILE *pczf, FILE *fout, int ip) {
  float *proj;

  proj = pczproj (pczf, ip);
  writeVector (proj, pczf->nframes, fout);
  sfree (proj);
}

/* This method dumps the average structure to a file in Amber format*/
void dumpavgx (PCZFILE *pczf, FILE *fout, int *atomList, int nAtoms) {
  float *avg, *selection;

  avg = pczavg (pczf);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packAtoms (avg, atomList, nAtoms);
    sfree (avg);
    avg = selection;
  }

  fprintf (fout, "%s\n", trim(pczf->title));
  writeVector (avg, nAtoms * 3, fout);
  sfree (avg);
}

/* This method dumps the atomic fluctuations to a file */
void dumpfluc (PCZFILE *pczf, FILE *fout, int iv, int bfactorsRequested, int *atomList, int nAtoms) {
  float *fluc, *selection;

  fluc = pczfluc (pczf, iv);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packValues (fluc, atomList, nAtoms);
    sfree (fluc);
    fluc = selection;
  }
  if (bfactorsRequested)
    pczbfactors (fluc, nAtoms);
  writeColumnVector (fluc, nAtoms, fout);
  sfree (fluc);
}

/* This method dumps the atomic fluctuations to a file in PDB format */
void dumpflucpdb (PCZFILE *pczf, FILE *fout, int iv, int bfactorsRequested, int *atomList, int nAtoms) {
  float *avg, *selection;
  float *fluc, *fselection;

  avg = pczavg (pczf);
  fluc = pczfluc (pczf, iv);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packAtoms (avg, atomList, nAtoms);
    fselection = packValues (fluc, atomList, nAtoms);
    sfree (avg);
    sfree (fluc);
    avg = selection;
    fluc = fselection;
  }

  if (bfactorsRequested)
    pczbfactors (fluc, nAtoms);
  dumppdb (pczf, fout, avg, fluc, atomList, nAtoms);
  sfree (avg);
  sfree (fluc);
}

/* This method dumps the RMS values between a given frame and all the others
 * to a file */
void dumprms (PCZFILE *pczf, FILE *fout, int ir) {
  int i, j;
  float *s1, *s2, rmsd;

  s2 = NULL;
  if (ir >= 0)
    s2 = pczscores (pczf, ir);

  for (i = 0; i < pczf->nframes; i++) {
    rmsd = 0.0;
    s1 = pczscores (pczf, i);
    if (ir >= 0)
      for (j = 0; j < pczf->nvecs; j++)
        s1[j] = s1[j] - s2[j];
    for (j = 0; j < pczf->nvecs; j++)
      rmsd += s1[j]*s1[j];
    rmsd /= pczf->natoms;
    fprintf (fout, "%8.3f\n", sqrt(rmsd));
    sfree (s1);
  }

  if (ir >= 0)
    sfree (s2);
}

/* This method dumps an animation along an eigenvector to a file
 * in Amber format */
void dumpanimx (PCZFILE *pczf, FILE *fout, int iv, int *atomList, int nAtoms) {
  float *avg, *evec, *proj, deform[20], x;
  float rmin, rmax, rinc;
  int i, j, colElems;
  float *selection, *eselection;

  avg = pczavg (pczf);
  evec = pczevec (pczf, iv);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packAtoms (avg, atomList, nAtoms);
    eselection = packAtoms (evec, atomList, nAtoms);
    sfree (avg);
    sfree (evec);
    avg = selection;
    evec = eselection;
  }

  proj = pczproj (pczf, iv);

  rmin = minval (proj, pczf->nframes);
  rmax = maxval (proj, pczf->nframes);
  rinc = (float)((rmax - rmin) * 0.1);
  deform[0] = 0.0;
  for (i = 1; i < 6; i++)
    deform[i] = deform[i - 1] + rinc;
  for (i = 6; i < 16; i++)
    deform[i] = deform[i - 1] - rinc;
  for (i = 16; i < 20; i++)
    deform[i] = deform[i - 1] + rinc;

  fprintf (fout, "Animation of eigenvector %i\n", iv + 1);
  for (j = 0; j < 20; j++) {
    colElems = 0;
    for (i = 0; i < nAtoms * 3; i++) {
      x = avg[i] + evec[i] * deform[j];
      if (colElems == 10) {
        colElems = 0;
        fprintf (fout, "\n");
      }
      fprintf (fout, "%8.3f", x);
      colElems++;
    }
    fprintf (fout, "\n");
  }

  sfree (proj);
  sfree (evec);
  sfree (avg);
}

/* This method dumps the average structure to a file in PDB format*/
void dumpavgpdb(PCZFILE *pczf, FILE *fout, int *atomList, int nAtoms/*, char **pdbline, int nr*/) {
  float *avg, *selection;

  avg = pczavg (pczf);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packAtoms (avg, atomList, nAtoms);
    sfree (avg);
    avg = selection;
  }

  dumppdb (pczf, fout, avg, NULL, atomList, nAtoms);
  sfree (avg);
}

/* This method dumps the eigenvectors to a file in a PDB-like formmatted file*/
void dumpevecpdb(PCZFILE *pczf, FILE *fout, int nEvec, int *atomList, int nAtoms) {
  float *evec, *selection;

  evec = pczevec (pczf, nEvec);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packAtoms (evec, atomList, nAtoms);
    sfree (evec);
    evec = selection;
  }

  dumppdb (pczf, fout, evec, NULL, atomList, nAtoms);
  sfree (evec);
}

/* This method dumps an animation along an eigenvector to a file
 * in PDB format */
void dumpanimpdb (PCZFILE *pczf, FILE *fout, int iv, int *atomList, int nAtoms/*, char **pdbline, int nr*/) {
  float *avg, *evec, *proj, deform[20], *snapshot;
  float rmin, rmax, rinc;
  int i, j;
  float *selection, *eselection;

  avg = pczavg (pczf);
  evec = pczevec (pczf, iv);
  if (atomList == NULL) {
    nAtoms = pczf->natoms;
  } else {
    selection = packAtoms (avg, atomList, nAtoms);
    eselection = packAtoms (evec, atomList, nAtoms);
    sfree (avg);
    sfree (evec);
    avg = selection;
    evec = eselection;
  }

  proj = pczproj (pczf, iv);
  snapshot = (float *)smalloc (pczf->natoms * 3 * sizeof (float));

  rmin = minval (proj, pczf->nframes);
  rmax = maxval (proj, pczf->nframes);
  rinc = (float)((rmax - rmin) * 0.1);
  deform[0] = 0.0;
  for (i = 1; i < 6; i++)
    deform[i] = deform[i - 1] + rinc;
  for (i = 6; i < 16; i++)
    deform[i] = deform[i - 1] - rinc;
  for (i = 16; i < 20; i++)
    deform[i] = deform[i - 1] + rinc;

  fprintf (fout, "REMARK Animation of eigenvector %4i\n", iv + 1);
  for (j = 0; j < 20; j++) {
    fprintf (fout, "MODEL  %5i\n", j+1);

    for (i = 0; i < nAtoms * 3; i++) {
      snapshot[i]   = avg[i]   + evec[i]   * deform[j];
    }
    dumppdb (pczf, fout, snapshot, NULL, atomList, nAtoms);

    fprintf (fout, "ENDMDL\n");
  }

  sfree (snapshot);
  sfree (proj);
  sfree (evec);
  sfree (avg);
}

/* Method to show the header information */
void dumpinfo (PCZFILE *pczf, FILE *fout) {
  fprintf (fout, "Title             : %s\n",      trim(pczf->title));
  fprintf (fout, "Atoms             : %8i\n",     pczf->natoms);
  fprintf (fout, "Vectors           : %8i\n",     pczf->nvecs);
  fprintf (fout, "Frames            : %8i\n",     pczf->nframes);
  fprintf (fout, "Total variance    : %8.2f\n",   pczf->vtot);
  fprintf (fout, "Explained variance: %8.2f\n",   pczf->vexplained);
  fprintf (fout, "Quality           : %8.2f%%\n", pczf->quality);
  fprintf (fout, "Dimensionality    : %8i\n",     pczf->dimensionality);
  fprintf (fout, "RMSd type         : ");
  switch (pczf->rms) {
  case STD_RMS:
    fprintf (fout, "Standard RMSd\n");
    break;
  case GAUSSIAN_RMS:
    fprintf (fout, "Gaussian RMSd\n");
    break;
  case NONE_RMS:
    fprintf (fout, "Unknown RMSd (if any)\n");
    break;
  default:
    fprintf (fout, "Unknown RMSd\n");
    break;
  }
  fprintf (fout, "Have atom names   : ");
  switch (pczf->haveAtomNames) {
  case FALSE:
    fprintf (fout, "False\n");
    break;
  default:
    fprintf (fout, "True\n");
    break;
  }
}

/* This method dumps the Lindemann coefficient for the selectetd atoms
 * to a file */
void dumplindemann (PCZFILE *pczf, FILE *fout, int *atomList, int nAtoms) {
  fprintf (fout, "%8.3f\n", pczlindemann (pczf, nAtoms, atomList));
}

/* This method dumps the force constants for each residue to a file */
void dumpforcecte (PCZFILE *pczf, FILE *fout, float temperature) {
  float *forcectes;

  forcectes = pczforcecte (pczf, temperature);
  writeColumnVector (forcectes, pczf->nvecs, fout);

  sfree (forcectes);
}

/* Method to safely compute runnning averages on the fly */
float runningAverage (float *vector, int size, int position, int howmany) {
  int i, nElems;
  float accumulator;
  
  accumulator = 0.0;
  nElems = 0;
  for (i = position-howmany; i <= position+howmany; i++) {
    if (i < 0 || i >= size)
      continue;
    else {
      accumulator += vector[i];
      nElems++;
    }
  }
  
  return accumulator/(float)nElems;
}

/* This method tries to predict the protein hinge points using atomic
 * fluctuation information */
void calchinge_fluc (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList) {
  float *fluc, *allFluc, **slope, maxFluc, threshold;
  float downVariance, upVariance, downMean, upMean;
  float *neighboursFluctuation, preMean, postMean, minFluct, maxFluct;
  int i, j, window, minWindow, maxWindow, checkingWindow;
  int coherent, downCounter, upCounter, downCandidates, upCandidates;
  int nNeighbours, hingePoint, preResidues, postResidues;
  int runAvgWide = 5;
  int firstResidueNumber = pczf->atomNames[atomList[0]].residueNumber;

  fprintf (fout, "#### Distance variation method\n");

  maxWindow = (int)roundf ((float)nResidues / 3.0);
  minWindow = maxWindow / 2;
  /*checkingWindow = minWindow;*/
  checkingWindow = 15;

  /* Now we extract the atomic fluctuations corresponding to the atoms
   * selected. I differentiate if I need the trajectory or only one
   * eigenvector */
  allFluc = pczfluc (pczf, nEvec);

  fluc = packValues (allFluc, atomList, nResidues);
  sfree (allFluc);

  /* We search the highest fluctuation in order to put a proper threshold */
  maxFluc = fluc[0];
  for (i = 1; i < nResidues; i++) {
    maxFluc = maxFluc<fabs(fluc[i]) ? fabs(fluc[i]) : maxFluc;
  }
  threshold = 0.05 * maxFluc;
  printf ("maxFluc=%f, threshold=%f\n", maxFluc, threshold);

  /* Calculus of the slope with different windows */
  slope = (float **)smalloc ((maxWindow - minWindow + 1) * sizeof (float *));
  for (window = minWindow; window <= maxWindow; window++) {
    slope[window-minWindow] = (float *)smalloc ((nResidues - window + 1) * sizeof (float));
    for (i = 0; i < nResidues - window + 1; i++) {
      /*slope[window-minWindow][i] = fluc[i] - fluc[i+window-1];*/
      slope[window-minWindow][i] = runningAverage (fluc, nResidues, i, runAvgWide) -
        runningAverage (fluc, nResidues, i+window-1, runAvgWide);
    }
  }

  /*
  for (window = minWindow; window <= maxWindow; window++) {
    printf ("#Window %i\n", window);
    for (i = 0; i < nResidues - window + 1; i++) {
      printf ("%5i %8.3f\n", i, slope[window-minWindow][i]);
    }
    printf ("\n");
  }
  //*/

  /* Search of the disruption point with the different windows */
  downCounter = 0;
  upCounter = 0;
  downCandidates = 0;
  upCandidates = 0;
  downVariance = 0.0;
  upVariance = 0.0;
  for (window = minWindow; window <= maxWindow; window++) {
    for (i = 0; i < nResidues - window + 1; i++) {
      /* Check for downhill hinge */
      if (slope[window-minWindow][i] > threshold) {
        coherent = TRUE;
        for (j = 1; j < checkingWindow; j++) {
          if (i + j >= nResidues - window + 1)
            continue;
          else
            coherent = coherent && (slope[window-minWindow][i+j] < threshold);
        }
        for (j = -checkingWindow; j < 0; j++) {
          if (i + j < 0)
            continue;
          else
            coherent = coherent && (slope[window-minWindow][i+j] > threshold);
        }
        if (coherent) {
          fprintf (fout, "Window %2i: residue %3i seems a downhill hinge point\n", window, pczf->atomNames[atomList[i]].residueNumber);
          downCounter++;
          downCandidates += i + firstResidueNumber;
          downVariance += (float)(i + firstResidueNumber)*(float)(i + firstResidueNumber);
        }
      }
      /* Check for uphill hinge */
      if (slope[window-minWindow][i] < -threshold) {
        coherent = TRUE;
        for (j = 1; j < checkingWindow; j++) {
          if (i + j >= nResidues - window + 1)
            continue;
          else
            coherent = coherent && (slope[window-minWindow][i+j] < -threshold);
        }
        for (j = -checkingWindow; j < 0; j++) {
          if (i + j < 0)
            continue;
          else
            coherent = coherent && (slope[window-minWindow][i+j] > -threshold);
        }
        if (coherent) {
          fprintf (fout, "Window %2i: residue %3i seems a uphill hinge point\n", window, pczf->atomNames[atomList[i + window - 1]].residueNumber);
          upCounter++;
          upCandidates += i + window + firstResidueNumber - 1;
          upVariance += (float)(i + window + firstResidueNumber - 1)*(float)(i + window + firstResidueNumber - 1);
        }
      }
    }
  }
  downMean     = (float)downCandidates/(float)downCounter;
  upMean       = (float)  upCandidates/(float)  upCounter;
  downVariance = (downVariance / (float)downCounter) - (downMean * downMean);
  upVariance   = (upVariance   / (float)upCounter)   - (upMean   * upMean);
  if (downCounter > 0) {
    fprintf (fout, "Consensus Downhill hinge point: %5.1f (%8.3f)\n", downMean, downVariance);
  }
  if (upCounter > 0) {
    fprintf (fout, "Consensus Uphill hinge point  : %5.1f (%8.3f)\n", upMean, upVariance);
  }

  /* Use of the neighbours distances to distinguish between a hinge point
   * and a unfold point */
  nNeighbours = 5;
  neighboursFluctuation = meanNeighboursDistanceFluctuations(pczf, nEvec, nResidues, atomList, nNeighbours);
  minFluct = neighboursFluctuation[0];
  maxFluct = neighboursFluctuation[0];
  for (i = 1; i < nResidues; i++) {
    if (neighboursFluctuation[i] < minFluct)
      minFluct = neighboursFluctuation[i];
    if (neighboursFluctuation[i] > maxFluct)
      maxFluct = neighboursFluctuation[i];
  }

  if (downCounter > 0) {
    /* The first residues are fluctuating more */
    hingePoint = (int)roundf (downMean);
    preMean = 0.0;
    postMean = 0.0;
    preResidues = 0;
    postResidues = 0;
    for (i = 0; i < hingePoint - 3; i++) {
      if (i < 0 || i >= nResidues)
        continue;
      preMean += neighboursFluctuation[i];
      preResidues++;
    }
    preMean /= preResidues;
    for (i = hingePoint + 3; i < nResidues; i++) {
      if (i < 0 || i >= nResidues)
        continue;
      postMean += neighboursFluctuation[i];
      postResidues++;
    }
    postMean /= postResidues;

    if (preMean > postMean - 0.4*(maxFluct-minFluct) &&
      preMean < postMean + 0.4*(maxFluct-minFluct)) {
      fprintf (fout, "Downhill hinge point seems a real hinge point\n");
    } else {
      fprintf (fout, "Downhill hinge point seems an unfolding point\n");
    }
  }
  if (upCounter > 0) {
    /* The last residues are fluctuating more */
    hingePoint = (int)roundf (downMean);
    preMean = 0.0;
    postMean = 0.0;
    preResidues = 0;
    postResidues = 0;
    for (i = 0; i < hingePoint - 3; i++) {
      if (i < 0 || i >= nResidues)
        continue;
      preMean += neighboursFluctuation[i];
      preResidues++;
    }
    preMean /= preResidues;
    for (i = hingePoint + 3; i < nResidues; i++) {
      if (i < 0 || i >= nResidues)
        continue;
      postMean += neighboursFluctuation[i];
      postResidues++;
    }
    postMean /= postResidues;

    if (preMean > postMean - 0.4*(maxFluct-minFluct) &&
      preMean < postMean + 0.4*(maxFluct-minFluct)) {
      fprintf (fout, "Uphill hinge point seems a real hinge point\n");
    } else {
      fprintf (fout, "Uphill hinge point seems an unfolding point\n");
    }
  }

  /* Resource freeing */
  sfree (neighboursFluctuation);
  for (window = minWindow; window <= maxWindow; window++) {
    sfree (slope[window-minWindow]);
  }
  sfree (slope);
  sfree (fluc);
}


/* This method tries to predict the protein hinge points checking the
 * distance variations between atoms */
void calchinge_distances (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList) {
  int i, j, indexPeak;

  float *avgDistances, *snapshot, *selection, *distances;
  float *sqdiff, diff, peak;

  fprintf (fout, "#### Force constant method\n");

  /* Now we can check distances */
  avgDistances = (float *)smalloc (nResidues * sizeof (float));
  distances    = (float *)smalloc (nResidues * sizeof (float));
  sqdiff       = (float *)smalloc (nResidues * sizeof (float));

  for (j = 0; j < nResidues; j++) {
    sqdiff[j] = 0.0;
    avgDistances[j] = 0.0;
  }

  /* First we compute the average distances */
  for (i = 0; i < pczf->nframes; i++) {
    snapshot = getFrameCoordinates (pczf, i, nEvec);
    selection = packAtoms (snapshot, atomList, nResidues);
    sfree (snapshot);
    getAverageDistances (selection, nResidues, distances);
    sfree (selection);
    for (j = 0; j < nResidues; j++) {
      avgDistances[j] += distances[j];
    }
  }
  for (j = 0; j < nResidues; j++) {
    avgDistances[j] /= pczf->nframes;
  }

  /* Here we compute the difference between average and frame distances */
  for (i = 0; i < pczf->nframes; i++) {
    snapshot = getFrameCoordinates (pczf, i, nEvec);
    selection = packAtoms (snapshot, atomList, nResidues);
    sfree (snapshot);
    getAverageDistances (selection, nResidues, distances);
    sfree (selection);
    for (j = 0; j < nResidues; j++) {
      diff = distances[j] - avgDistances[j];
      sqdiff[j] += (diff * diff);
    }
  }

  for (j = 0; j < nResidues; j++) {
    sqdiff[j] = 1.0 / (sqdiff[j] / pczf->nframes);
  }

  /* Search for a peak */
  indexPeak = 0;
  peak = sqdiff[0];
  for (j = 1; j < nResidues; j++) {
    if (peak < sqdiff[j]) {
      indexPeak = j;
      peak = sqdiff[j];
    }
  }

  /* Finally print the results */
  fprintf (fout, "# Peak constant (possible hinge): residue %i (%.3f)\n", pczf->atomNames[atomList[indexPeak]].residueNumber, peak);
  fprintf (fout, "# Force constants per residue\n");
  writeColumnVector (sqdiff, nResidues, fout);

  sfree (sqdiff);
  sfree (distances);
  sfree (avgDistances);
}

/* This method computed the apparent stiffness between residues */
void dumpstiffness (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList, float temperature) {
  int i, j, k;

  float **avgDistances, *snapshot, *selection, **distances;
  float **sqdiff, diff;
  float const kB = 0.00198717; /*1.3806503e-23*/

  if (atomList == NULL) {
    atomList = selectAtoms (pczf, ":1|~:1", &nResidues);
  }


  /* Now we can check distances */
  avgDistances = (float **)smalloc (nResidues * sizeof (float *));
  distances    = (float **)smalloc (nResidues * sizeof (float *));
  sqdiff       = (float **)smalloc (nResidues * sizeof (float *));

  for (j = 0; j < nResidues; j++) {
    distances[j] = (float *)smalloc (nResidues * sizeof (float));
    sqdiff[j] = (float *)smalloc (nResidues * sizeof (float));
    avgDistances[j] = (float *)smalloc (nResidues * sizeof (float));
    for (i = 0; i < nResidues; i++) {
      sqdiff[j][i] = 0.0;
      avgDistances[j][i] = 0.0;
      distances[j][i] = 0.0;
    }
  }

  /* First we compute the average distances */
  for (i = 0; i < pczf->nframes; i++) {
    snapshot = getFrameCoordinates (pczf, i, nEvec);
    selection = packAtoms (snapshot, atomList, nResidues);
    sfree (snapshot);
    for (j = 0; j < nResidues; j++) {
      for (k = j; k < nResidues; k++) {
        avgDistances[j][k] += distance (selection, j, k);
      }
    }
    sfree (selection);
  }
  for (j = 0; j < nResidues; j++) {
    for (k = j; k < nResidues; k++) {
      avgDistances[j][k] /= pczf->nframes;
    }
  }

  /* Here we compute the difference between average and frame distances */
  for (i = 0; i < pczf->nframes; i++) {
    snapshot = getFrameCoordinates (pczf, i, nEvec);
    selection = packAtoms (snapshot, atomList, nResidues);
    sfree (snapshot);
    for (j = 0; j < nResidues; j++) {
      for (k = j; k < nResidues; k++) {
        diff = distance (selection, j, k) - avgDistances[j][k];
        sqdiff[j][k] += (diff * diff);
      }
    }
    sfree (selection);
  }

  for (j = 0; j < nResidues; j++) {
    for (k = j; k < nResidues; k++) {
      sqdiff[j][k] = (kB * temperature) / (sqdiff[j][k] / pczf->nframes);
      sqdiff[k][j] = sqdiff[j][k];
    }
  }
  for (j = 0; j < nResidues; j++) {
    sqdiff[j][j] = 0.0;
  }
  
  /* Finally print the results */
  for (i = 0; i < nResidues; i++) {
    writeCSVSci (sqdiff[i], nResidues, fout);
  }

  for (j = 0; j < nResidues; j++) {
    sfree (distances[j]);
    sfree (sqdiff[j]);
    sfree (avgDistances[j]);
  }
  sfree (sqdiff);
  sfree (distances);
  sfree (avgDistances);
}

/* This method computes the distance between two points.
 * i and j indexes start at zero */
float distance (float *xyz, int i, int j) {
	float xi, yi, zi, xj, yj, zj, dx, dy, dz, dist;
	
	xi = xyz[i*3];
	yi = xyz[i*3+1];
	zi = xyz[i*3+2];
	xj = xyz[j*3];
	yj = xyz[j*3+1];
	zj = xyz[j*3+2];
	
	dx = xi - xj;
	dy = yi - yj;
	dz = zi - zj;
	
	dist = sqrt (dx*dx + dy*dy + dz*dz);
	
	return dist;
}

float clusterDistance (ClusterList *cl, int clusterA, int clusterB, float **drij) {
  int elemA, elemB, i, j;
  float max;
  
  max = 0.0;
  for (i = 0; i < getNElemsOfCluster (cl, clusterA); i++) {
    elemA = getElemOfCluster (cl, clusterA, i);
    for (j = 0; j < getNElemsOfCluster (cl, clusterB); j++) {
      elemB = getElemOfCluster (cl, clusterB, j);
      max = max>drij[elemA][elemB] ? max : drij[elemA][elemB];
    }
  }
  
  return max;
}

float hierarchicalClusterDistance (ClusterList *cl, int clusterA, int clusterB, float **drij) {
  int elemA, elemB;
  float max;
  IntegerList *ilA, *ilB, *ptrA, *ptrB;
  
  ilA = getListForCluster (cl, clusterA);
  ilB = getListForCluster (cl, clusterB);
  
  max = 0.0;
  for (ptrA = ilA; ptrA != NULL; ptrA = ptrA->next) {
    elemA = ptrA->value;
    for (ptrB = ilB; ptrB != NULL; ptrB = ptrB->next) {
      elemB = ptrB->value;
      max = max>drij[elemA][elemB] ? max : drij[elemA][elemB];
    }
  }
  
  deleteIntegerList (ilA);
  deleteIntegerList (ilB);
  
  return max;
}

int allBelowThreshold (ClusterList *cl, int clusterId, int elem, float **drij, float threshold) {
  int below, i, elemCluster;
  
  below = TRUE;
  for (i = 0; below && (i < getNElemsOfCluster (cl, clusterId)); i++) {
    elemCluster = getElemOfCluster (cl, clusterId, i);
    below = below &&
      drij[elem][elemCluster] < threshold &&
      drij[elemCluster][elem] < threshold;
    /*
    if (!below) {
      fprintf (stdout, "Dos elements sobre el threshold: %i i %i (%f)\n", elem, elemCluster, drij[elem][elemCluster]);
    }
    */
  }
  
  return below;
}

float meanDistanceElemToCluster (ClusterList *cl, int clusterId, float **drij, int elem) {
  float distance;
  int nElems, i;

  distance = 0.0;
  nElems = 0;  
  for (i = 0; i < getNElemsOfCluster (cl, clusterId); i++) {
    if (elem != getElemOfCluster (cl, clusterId, i)) {
      distance += drij[elem][getElemOfCluster (cl, clusterId, i)];
      nElems++;
    }
  }
  
  distance /= (float)nElems;
  
  return distance;
}

/* This method tries to predict the protein hinge points checking the
 * domains predicted with the method of Lavery's paper */
void calchinge_domains (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList) {
  int i, j, k, mini, minj, trobat, clusterFound, clusterId, testCluster, elemId;
  float *allmean, *mean, **drij, *snapshot, *selection;
  float min, lastmin, meanA, meanB, cdist, minimumDistanceBtwClusters;
  float mindrij, maxdrij;
  /*float threshold;*/
  int elem, migrationPerformed, clustersJoined, nIterations;
  ClusterList *cl;
  int clid, clid2, atid, atid2, *interactions, *hinges, idxInter, idxHinge;
  double interactionDistance;
  int firstResidueNumber = pczf->atomNames[atomList[0]].residueNumber;
  
  fprintf (fout, "#### Lavery method\n");
  
  allmean = pczavg (pczf);
  mean = packAtoms (allmean, atomList, nResidues);
  sfree (allmean);
  
  /* Initialization of the delta-r(i,j) matrix */
  drij = (float **)smalloc (nResidues * sizeof (*drij));
  for (i = 0; i < nResidues; i++) {
    drij[i] = (float *)smalloc (nResidues * sizeof (**drij));
    for (j = 0; j < nResidues; j++) {
      drij[i][j] = 0.0;
    }
  }
  
  fprintf (fout, "*** Evec used: %i\n", nEvec);
  fprintf (fout, "*** Threshold defined: %f\n", threshold);
  /* We now compute the average distance oscillations */
  for (k = 0; k < pczf->nframes; k++) {
    snapshot = getFrameCoordinates (pczf, k, nEvec);
    selection = packAtoms (snapshot, atomList, nResidues);
    sfree (snapshot);

    for (i = 0; i < nResidues; i++) {
      for (j = 0; j < nResidues; j++) {
        drij[i][j] += fabs (distance (mean, i, j) - distance (selection, i, j));
      }
    }
    sfree (selection);
  }
  mindrij = drij[1][0]/(float)pczf->nframes;
  maxdrij = drij[1][0]/(float)pczf->nframes;
  for (i = 0; i < nResidues; i++) {
    for (j = 0; j < nResidues; j++) {
      drij[i][j] /= (float)pczf->nframes;
      mindrij = (mindrij>drij[i][j]&&drij[i][j]>0.0)?drij[i][j]:mindrij;
      maxdrij = maxdrij<drij[i][j]?drij[i][j]:maxdrij;
    }
  }
  fprintf (fout, "*** Min. drij: %f\n", mindrij);
  fprintf (fout, "*** Max. drij: %f\n", maxdrij);

  /* Begin clustering process */
  
  cl = NULL;
  while (getNClusters(cl) >= 5 || getNClusters(cl) <= 1) {
    threshold = (mindrij + maxdrij) / 2.0;
    fprintf (stderr, "*** Trying threshold %f\n", threshold);
    
  
    
    
    
    
    clusterFound = TRUE;
    /*threshold = 0.2;*/
    cl = newClusterList (nResidues);
    while (clusterFound) {
      clusterFound = FALSE;
      /* Look for the next minimum not in any cluster */
      mini = -1;
      minj = -1;
      min = threshold;
      for (i = 0; i < nResidues; i++) {
        if (isElemInClusterList (cl, i))
          continue; /* Discard clusterized elements */
        for (j = 0; j < nResidues; j++) {
          if (isElemInClusterList (cl, j))
            continue; /* Discard clusterized elements */
          if (drij[i][j] > 0.0 && drij[i][j] < min) {
            mini = i;
            minj = j;
            min = drij[mini][minj];
          }
        }
      }
      if (mini != -1 && minj != -1) {
        /* Now generate the cluster for this domain */
        /*fprintf (fout, "El minim esta a i:%i j:%i i val %f\n", mini, minj, min);*/
        clusterFound = TRUE;
        clusterId = newCluster (cl);
        addElemToCluster (cl, clusterId, mini);
        addElemToCluster (cl, clusterId, minj);
        trobat = TRUE;
        while (trobat) {
          trobat = FALSE;
          lastmin = min;
          min = threshold;
          for (j = 0; j < nResidues; j++) {
            if (isElemInClusterList (cl, j))
              continue; /* Discard clusterized elements */
            if (drij[mini][j] > lastmin && drij[mini][j] < min) {
              minj = j;
              min = drij[mini][minj];
              trobat = TRUE;
            }
          }
          if (trobat) {
            /* We check that the new element is below the threshold for all the
             * existent elements in the cluster */
            if (allBelowThreshold (cl, clusterId, minj, drij, threshold)) {
              addElemToCluster (cl, clusterId, minj);
              /*fprintf (fout, "Afegim al cluster el residu %i (%f)\n", minj, drij[mini][minj]);*/
            } else {
              /*fprintf (fout, "NO Afegim al cluster el residu %i (%f)\n", minj, drij[mini][minj]);*/
            }
          }
        }
      }
    }
  
    /*printClusterList (cl, fout);*/
  
    /* Primary clusters generated. Now we will check that each residue is in the
     * cluster with the lowest mean distance with itself */
    if (getNClusters (cl) > 0) {
      migrationPerformed = TRUE;
      nIterations = 10;
      while (migrationPerformed && nIterations > 0) {
        /*fprintf (fout, "Looking for migrating elements\n");*/
        migrationPerformed = FALSE;
        nIterations--;
        for (clusterId = 0; clusterId < getNClusters (cl); clusterId++) {
          /*fprintf (fout, "Checking cluster %i\n", clusterId);*/
          for (elemId = 0; elemId < getNElemsOfCluster (cl, clusterId); elemId++) {
            /*fprintf (fout, "  Checking element %i\n", elemId);*/
            elem = getElemOfCluster(cl, clusterId, elemId);
            for (testCluster = 0; testCluster < getNClusters (cl); testCluster++) {
              /*fprintf (fout, "    Against cluster %i\n", testCluster);*/
              if (testCluster == clusterId)
                continue;
              /* We check if the element is below the threshold for the target cluster */
              if (allBelowThreshold (cl, testCluster, elem, drij, threshold)) {
                /*fprintf (fout, "    The element %i(cl:%i) could be migrated to cl:%i\n", elem, clusterId, testCluster);*/
                /* Now we compute the mean distance to decide if we must migrate the element */
                meanA = meanDistanceElemToCluster (cl, clusterId, drij, elem);
                meanB = meanDistanceElemToCluster (cl, testCluster, drij, elem);
                if (meanB < meanA) {
                  /*fprintf (fout, "      The change of cluster must be accepted %i (cl:%i->%i)(%f->%f)\n", elem, clusterId, testCluster, meanA, meanB);*/
                  deleteElemFromCluster (cl, clusterId, elem);
                  addElemToCluster (cl, testCluster, elem);
                  migrationPerformed = TRUE;
                  break; /* We just moved it, we must not check for another movement */
                }
              }
            }
          }
        }
      }
      if (nIterations <= 0) {
        fprintf (stderr, "Warning: too many iterations while refining domains\n");
      }
    
      /*printClusterList (cl, fout);*/
    
      /* We now join too near clusters */
      /*fprintf (fout, "--- Joining nearby clusters --- \n");*/
      minimumDistanceBtwClusters = threshold * 1.3; 
      clustersJoined = TRUE;
      while (clustersJoined) {
        clustersJoined = FALSE;
        newHierarchicalCluster (cl);
        for (i = 0; i < getNElemsInTopClusters (cl) && clustersJoined == FALSE; i++) {
          for (j = i; j < getNElemsInTopClusters (cl) && clustersJoined == FALSE; j++) {
            if (i == j)
              continue;
            cdist = hierarchicalClusterDistance (cl, getElemInTopClusters (cl, i), getElemInTopClusters (cl, j), drij);
            /*fprintf (fout, "Distance between clusters %i-%i -> %f\n", getElemInTopClusters (cl, i), getElemInTopClusters (cl, j), cdist);*/
            if (cdist < minimumDistanceBtwClusters) {
              /*fprintf (fout, "  ** Aquests dos clusters estan massa junts\n");*/
              joinClusters (cl, getElemInTopClusters (cl, i), getElemInTopClusters (cl, j));
              /*printClusterList (cl, fout);*/
              clustersJoined = TRUE;
            }
          }
        }
      }
      
      fprintf (fout, "nClusters: %i threshold: %f\n", getNClusters (cl), threshold);
      
      /* Now we are goint to build a hierarchical cluster tree */
      /*fprintf (fout, "--- Building hierarchy --- \n");*/
      while (getNElemsInTopClusters (cl) > 1) {
        mini = -1;
        minj = -1;
        for (i = 0; i < getNElemsInTopClusters (cl); i++) {
          for (j = i; j < getNElemsInTopClusters (cl); j++) {
            if (i == j)
              continue;
            cdist = hierarchicalClusterDistance (cl, getElemInTopClusters (cl, i), getElemInTopClusters (cl, j), drij);
            /*fprintf (fout, "Distance between clusters %i-%i -> %f\n", getElemInTopClusters (cl, i), getElemInTopClusters (cl, j), cdist);*/
            if (mini == -1 || cdist < min) {
                mini = i;
                minj = j;
                min = cdist;
            }
          }
        }
        /*fprintf (fout, "Minimum distance %f (%i, %i)\n", min, getElemInTopClusters (cl, mini), getElemInTopClusters (cl, minj));*/
        addHierarchicalCluster (cl, getElemInTopClusters (cl, mini), getElemInTopClusters (cl, minj));
      }
    }

    
    
    
    if (getNClusters(cl) >= 5 || getNClusters(cl) <= 1) {
      if (getNClusters(cl) >= 5) {
        mindrij = threshold;
      } else if (getNClusters(cl) <= 1) {
        maxdrij = threshold;
      }
      deleteClusterList (cl);
      cl = NULL;
    }
  }
    
    
    
    
    
    
    
    
    
  printClusterList (cl, fout, firstResidueNumber);

  /* Now we are going to find the interaction regions and the hinges */
  /* Tho C-alpha will be considered interacting if they are closer than 4
   * angstroms and pertain to different clusters */
  if (getNClusters(cl) == 1) {
    fprintf (fout, "# Only one cluster found. No hinge possible.\n");
  } else {
    interactionDistance = 4.0;
    interactions = (int *)smalloc (nResidues * sizeof (*interactions));
    idxInter = 0;
    for (clid = 0; clid < getNClusters(cl) - 1; clid++) {
      for (atid = 0; atid < getNElemsOfCluster(cl, clid); atid++) {
        for (clid2 = clid + 1; clid2 < getNClusters(cl); clid2++) {
          for (atid2 = 0; atid2 < getNElemsOfCluster(cl, clid2); atid2++) {
            if (distance(mean, getElemOfCluster (cl, clid, atid),
                getElemOfCluster (cl, clid2, atid2)) < interactionDistance) {
              idxInter = addToVector (interactions, getElemOfCluster (cl, clid, atid), idxInter);
              idxInter = addToVector (interactions, getElemOfCluster (cl, clid2, atid2), idxInter);
            }
          }
        }
      }
    }
    fprintf (fout, "Interacting regions:");
    for (i = 0; i < idxInter; i++) {
      fprintf (fout, " %i", interactions[i] + firstResidueNumber);
    }
    fprintf (fout, "\n");
    
    /* Now we look for the hinges: pairs of residues one next to another
     * found in the interacting regions */
    hinges = (int *)smalloc (nResidues * sizeof (*interactions));
    idxHinge = 0;
    for (i = 0; i < idxInter - 1; i++) {
      if (interactions[i] == interactions[i + 1] - 1) {
        idxHinge = addToVector (hinges, interactions[i], idxHinge);
        idxHinge = addToVector (hinges, interactions[i + 1], idxHinge);
      }
    }
    fprintf (fout, "Hinge residues:");
    for (i = 0; i < idxHinge; i++) {
      fprintf (fout, " %i", hinges[i] + firstResidueNumber);
    }
    fprintf (fout, "\n");
    
    sfree (hinges);
    sfree (interactions);
  }

  /* Resource freeing */
  deleteClusterList (cl);
  sfree (mean);
  for (i = 0; i < nResidues; i++) {
    sfree (drij[i]);
  }
  sfree (drij);
  
}

/* This method dumps the information given by the hinge prediction methods */
void dumphinge (PCZFILE *pczf, FILE *fout, int nEvec) {
  int *atomList, nResidues;

  /* First we check that we have all the information necessary to do the
   * calculus */
  if (!pczf->haveAtomNames) {
    fprintf (stderr, "Error: No atom names are present\n");
    return;
  }
  if (pczf->rms != GAUSSIAN_RMS) {
    fprintf (stderr, "Warning: Compression done without gaussian RMSd\n");
  }

  /* We select Calpha atoms only */
  atomList = selectAtoms (pczf, "@CA", &nResidues);

  if (nResidues == 0) {
    fprintf (stderr, "Error: There are no alpha carbons\n");
    return;
  } else {
    printf ("There are %i residues\n", nResidues);
  }

  /* Now we call the real computing routines with the CA atom list */
  /*****************************************************************/
  calchinge_fluc (pczf, fout, nEvec, nResidues, atomList);
  calchinge_distances (pczf, fout, nEvec, nResidues, atomList);
  calchinge_domains (pczf, fout, nEvec, nResidues, atomList);
  /*sfree (meanNeighboursDistanceFluctuations(pczf, nEvec, nResidues, atomList, 5));*/

  sfree (atomList);
}

/* This method computes the mean fluctuation of the neighbour's distances */
float *meanNeighboursDistanceFluctuations (PCZFILE *pczf, int nEvec, int nResidues, int *atomList, int nNeighbours) {
  float *avg, *allAvg, *allXyz, *xyz, **means, *meanOfMeans, *distanceFluctuation;
  float aux;
  int **neighbours;
  int i, j;
  FILE *fout = stderr;

  /* Now we extract the average corresponding to the atoms selected. */
  allAvg = pczavg (pczf);

  avg = packAtoms (allAvg, atomList, nResidues);
  sfree (allAvg);

  /* We compute the neighbours from the average structure */
  neighbours = getNeighbours (avg, nResidues, nNeighbours);
  /* Show the neighbour list */
  /*
  for (j = 0; j < nResidues; j++) {
    printf ("Residu %3i: ", j);
    for (i = 0; i < nNeighbours; i++) {
      printf ("%3i, ", neighbours[j][i]);
    }
    printf ("\n");
  }
  */

  /* Now we compute the neighbour distance fluctuation */
  /* I must create a method to get the distances from one residue to the ones
   * in the neighbour list.
   * Then I must take a mean distance for each residue with all the neighbour's
   * distances.
   * Once we have the distance's means for each frame, we compute a global mean
   * for all the frames.
   * After all, we only have to get the fluctuation around the mean:
   * sqrt(sum((mean(fi)-mean(f))^2))
   */
  means = (float **)smalloc (nResidues * sizeof (float *));
  for (i = 0; i < nResidues; i++) {
    means[i] = (float *)smalloc (pczf->nframes * sizeof (float));
  }
  for (j = 0; j < pczf->nframes; j++) {
    allXyz = getFrameCoordinates (pczf, j, nEvec);
    xyz = packAtoms (allXyz, atomList, nResidues);
    sfree (allXyz);
    for (i = 0; i < nResidues; i++) {
      means[i][j] = getMeanDistance (xyz, i, neighbours[i], nNeighbours);
    }
    sfree (xyz);
  }
  for (i = 0; i < nResidues; i++) {
    sfree(neighbours[i]);
  }
  sfree (neighbours);

  /* We compute the mean of the means */
  meanOfMeans = (float *)smalloc (nResidues * sizeof (float));
  for (i = 0; i < nResidues; i++) {
    meanOfMeans[i] = 0.0;
    for (j = 0; j < pczf->nframes; j++) {
      meanOfMeans[i] += means[i][j];
    }
    meanOfMeans[i] /= pczf->nframes;
  }

  /* Finally we compute the distance fluctuation */
  distanceFluctuation = (float *)smalloc (nResidues * sizeof (float));
  for (i = 0; i < nResidues; i++) {
    distanceFluctuation[i] = 0.0;
    for (j = 0; j < pczf->nframes; j++) {
      aux = means[i][j] - meanOfMeans[i];
      distanceFluctuation[i] += aux * aux;
    }
    distanceFluctuation[i] = sqrt(distanceFluctuation[i]);
  }

  /* Debug printing */
  /*
  fprintf (fout, "Res. Mean.1 Mean.2 MMean Fluct\n");
  for (i = 0; i < nResidues; i++) {
    fprintf (fout, "%4i %6.2f %6.2f %5.2f %5.2f\n", i+1, means[i][0], means[i][1], meanOfMeans[i], distanceFluctuation[i]);
  }
  */

  /* Neighbours fluctuations printing */
  fprintf (fout, "# Neighbour distance fluctuation, per residue\n");
  for (i = 0; i < nResidues; i++) {
    fprintf (fout, "%3i %8.3f\n", i+pczf->atomNames[atomList[0]].residueNumber, distanceFluctuation[i]);
  }

  sfree (meanOfMeans);
  for (i = 0; i < nResidues; i++) {
    sfree(means[i]);
  }
  sfree (means);
  sfree (avg);

  return (distanceFluctuation);
}

/* This method returns the mean distance from a residue to its neighbours */
float getMeanDistance (float *xyz, int residue, int *neighbours, int nNeighbours) {
  int i;
  float mean;

  mean = 0.0;

  for (i = 0; i < nNeighbours; i++) {
    mean += getDistance (xyz, residue, neighbours[i]);
  }
  mean /= (float)nNeighbours;

  return mean;
}

/* This method computes the nearest neighbours to each residue */
int **getNeighbours (float *xyz, int nResidues, int nNeighbours) {
  struct dist_t **distances;
  int **neighbours, i, j;

  /* We get the neighbours list */
  distances = getNeighboursExtended (xyz, nResidues, nNeighbours);

  /* Now we strip the distances from the neighbours list */
  neighbours = (int **)smalloc (nResidues * sizeof (int *));
  for (i = 0; i < nResidues; i++) {
    neighbours[i] = (int *)smalloc (nNeighbours * sizeof (int));
    for (j = 0; j < nNeighbours; j++) {
      neighbours[i][j] = distances[i][j].residue;
    }
  }

  /* Resource freeing */
  for (i = 0; i < nResidues; i++) {
    sfree (distances[i]);
  }
  sfree (distances);

  return neighbours;
}

/* This method gets the distances from each residue to its neighbours */
struct dist_t **getNeighboursExtended (float *xyz, int nResidues, int nNeighbours) {
  struct dist_t **distances;
  int i;

  distances = (struct dist_t **)smalloc (nResidues * sizeof (struct dist_t*));
  for (i = 0; i < nResidues; i++) {
    distances[i] = getDistances (xyz, nResidues, nNeighbours, i);
  }

  return distances;
}

/* This method gets the distance from one residue to its neighbours */
struct dist_t *getDistances (float *xyz, int nResidues, int nNeighbours, int residue) {
  struct dist_t *distances;
  int i;
  float distance;

  distances = (struct dist_t *)smalloc (nNeighbours * sizeof (struct dist_t));
  for (i = 0; i < nNeighbours; i++) {
    distances[i].residue = -1;
    distances[i].distance = -1.0;
  }

  for (i = 0; i < nResidues; i++) {
    if ((i == residue) || (i == residue - 1) || (i == residue + 1))
      continue;
    distance = getDistance (xyz, residue, i);
    addDistanceList (distances, nNeighbours, i, distance);
  }

  return distances;
}

/* This method inserts a new element into the ordered vector. The vector
 * contains a list of the nearest residues and their distances to the residue
 * which is associated to this list.
 * The ordered structure is called distances.
 * The number of neighbours that we are maintaining is nNeighbours
 * The residue number that we are introducing is residue
 * The distance associated to this residue is distance */
void addDistanceList (struct dist_t *distances, int nNeighbours, int residue, float distance) {
  int i, lastElement, position;

  if (distances[0].residue == -1) {
    /* initialization case: empty array */
    distances[0].residue = residue;
    distances[0].distance = distance;
  } else {
    /* Iterative case: at least one element */
    /* First we look for the last element */
    lastElement = nNeighbours - 1;
    while (distances[lastElement].residue == -1)
      lastElement--;
    /* Then we look for the position */
    position = lastElement;
    while (position >= 0 && distances[position].distance > distance)
      position--;
    if (position == (nNeighbours - 1)) {
      if (distances[position].distance > distance) {
        distances[position].residue = residue;
        distances[position].distance = distance;
      }
    } else {
      if (position >= 0) {
        /* Move the last elements and insert new element here */
        for (i = nNeighbours - 1; i > position + 1; i--) {
          distances[i].residue = distances[i-1].residue;
          distances[i].distance = distances[i-1].distance;
        }
        distances[position+1].residue = residue;
        distances[position+1].distance = distance;
      } else {
        /* Move all elements and insert new at beginning */
        for (i = nNeighbours - 1; i > 0; i--) {
          distances[i].residue = distances[i-1].residue;
          distances[i].distance = distances[i-1].distance;
        }
        distances[0].residue = residue;
        distances[0].distance = distance;
      }
    }
  }
}

/* Get the distance between two residues, residue1 and residue2 */
float getDistance (float *xyz, int residue1, int residue2) {
  float distance;
  int x1, y1, z1, x2, y2, z2;
  float dx, dy, dz;

  x1 = (residue1*3);
  y1 = (residue1*3)+1;
  z1 = (residue1*3)+2;
  x2 = (residue2*3);
  y2 = (residue2*3)+1;
  z2 = (residue2*3)+2;
  dx = xyz[x2] - xyz[x1];
  dy = xyz[y2] - xyz[y1];
  dz = xyz[z2] - xyz[z1];
  distance = sqrt(dx*dx + dy*dy + dz*dz);

  return distance;
}

/*
void getDistances (PCZFILE *pczf, int residue, int *atomList, int nResidues, struct dist_t *distances) {
  for (i = 0; i < pczf->nframes; i++) {
    snapshot = pczsnap (pczf, i);
    if (snapshot == NULL)
      return -1;
    writeVector (snapshot, pczf->natoms * 3, fout);
    sfree (snapshot);

}

int **getNeighbourList (PCZFILE *pczf, int *atomList, int nResidues) {
  int **neighbours;
  struct dist_t **distances;
  int i, nNeighbours;

  nNeighbours = 10;

  distances = (struct dist_t **)smalloc (nResidues * sizeof (struct dist_t*));
  for (i = 0; i < nResidues; i++) {
    distances[i] = (struct dist_t *)smalloc (nNeighbours * sizeof (struct dist_t));
    getNeighbours (pczf, i, atomList, nResidues, distances, nNeighbours);
  }

  return neighbours;
}

int compareDistances (const void *m1, const void *m2) {
  struct dist_t *d1 = (struct dist_t *)m1;
  struct dist_t *d2 = (struct dist_t *)m2;

  if (d1->distance < d2->distance)
    return -1;
  else if (d1->distance < d2->distance)
    return 1;
  else
    return 0;
}

*/

/* Given a snapshot, it computes the average distance between each residue
 * against all other residues (excluding its immediate neigbours, whose
 * distances are virtually maintained) */
void getAverageDistances (float *xyz, int nResidues, float *distances) {
  int i, j, consideredAtoms;
  float dx, dy, dz;

  for (i = 0; i < nResidues; i++) {
    distances[i] = 0.0;
    consideredAtoms = 0;
    for (j = 0; j < nResidues; j++) {
      if (j == i || j == i-1 || j == i+1)
        continue;
      dx = xyz[i*3  ] - xyz[j*3  ];
      dy = xyz[i*3+1] - xyz[j*3+1];
      dz = xyz[i*3+2] - xyz[j*3+2];
      distances[i] += sqrt (dx * dx + dy * dy + dz * dz);
      /*printf ("(%i, %i) dx: %8.3f, dy: %8.3f, dz: %8.3f, distance: %8.3f\n", i, j, dx, dy, dz, distances[i]);*/
      consideredAtoms++;
    }
    distances[i] /= consideredAtoms;
  }
}

/* This method computes the motion collectivity index for a given eigenvector */
void dumpcollectivity (PCZFILE *pczf, FILE *fout, int nEvecCollectivity) {
  float *evec, alpha, sum, normalizedValue, collectivity;
  float *u;
  int i, j, begin, end;

  if (nEvecCollectivity == -1) {
    /* Compute data for all eigenvectors */
	  begin = 0;
	  end   = pczf->nvecs - 1;
  } else {
    /* Compute data for only one eigenvector */
	  begin = nEvecCollectivity;
	  end   = nEvecCollectivity;
  }
  
  u = (float *)smalloc (pczf->natoms * sizeof (float));
  for (j = begin; j <= end; j++) {
	  alpha = 0.0;
	  evec = pczevec (pczf, j);
	  for (i = 0; i < pczf->natoms; i++) {
	    u[i] = evec[i*3] * evec[i*3] +
	           evec[i*3+1] * evec[i*3+1] +
	           evec[i*3+2] * evec[i*3+2];
	    alpha += u[i];
	  }
	  alpha = 1.0 / alpha;
	  sum = 0.0;
	  for (i = 0; i < pczf->natoms; i++) {
	    normalizedValue = u[i] * alpha;
	    sum += normalizedValue * log (normalizedValue);
	  }
	  collectivity = expf (-sum) / pczf->natoms;
	
	  fprintf (fout, "%f\n", collectivity);
	
	  sfree (evec);
  }
  sfree (u);
}

/* This method computes the Mahalanobis distance */
void dumpmahalanobis (PCZFILE *pczf, FILE *fout, int nev, char *xfile, int *atomList, int nAtoms) {
  pczmahalanobis (pczf, fout, nev, xfile, atomList, nAtoms);
}

/* Information method fot the user */
void helptext (void) {
  printf("Usage:\npczdump -i infile [-o outfile]\n");
  printf("\t--info           : Show general information about the file\n");
  printf("\t--avg            : Show the average structure\n");
  printf("\t--evals          : Show the eigenvalues\n");
  printf("\t--evec iv        : Show the <iv> eigenvector\n");
  printf("\t--proj iv        : Show the projections for <iv> eigenvector\n");
  printf("\t--rms iref       : Compute the RMS of the traj. against frame <iref>\n");
  printf("\t--fluc iv        : Compute atomic fluctuations for <iv> evec (0 for all)\n");
  printf("\t--bfactor        : Give the fluctuation values as B-factors\n");
  printf("\t--anim iv        : Animate the average structure along <iv> eigenvector\n");
  printf("\t--lindemann      : Compute Lindemann coefficient\n");
  printf("\t--collectivity iv: Compute a collectivity coefficient for <iv> evec\n");
  printf("\t--forcecte temp  : Compute force constants for temperature <temp>\n");
  printf("\t--hinge          : Look for hinge points (need a file with atom names\n");
  printf("\t                   and gaussian RMS applied\n");
  printf("\t--mahref         : Trajectory file to use in the Mahalanobis distance\n");
  printf("\t                   calculus\n");
  printf("\t--mahnev         : Number of eigenvectors to use in the Mahalanobis\n");
  printf("\t                   distance calculus\n");
  printf("\t--mask <mask>    : ptraj-like selection mask\n");
  printf("\t--pdb            : Ask for output in PDB format (if suitable)\n");
  printf("\t--verbose        : Show info about progress\n");
  printf("\t--help           : Show this message\n");
}

CLParameters *readCLParameters (int argc, char **argv) {
  static struct option long_options[] = {
    {"infile",        1, 0, 'i'},
    {"outfile",       1, 0, 'o'},
    {"info",          0, 0, 'I'},
    {"avg",           0, 0, 'a'},
    {"evec",          1, 0, 'e'},
    {"proj",          1, 0, 'p'},
    {"rms",           1, 0, 'r'},
    {"bfactor",       0, 0, 'b'},
    {"fluc",          2, 0, 'f'},
    {"evals",         0, 0, 'E'},
    {"anim",          1, 0, 'A'},
    {"pdb",           0, 0, 'P'},
    {"lindemann",     0, 0, 'l'},
    {"forcecte",      1, 0, 'c'},
    {"temperature",   1, 0, 'T'},
    {"collectivity",  2, 0, 'C'},
    {"hinge",         2, 0, 'H'},
    {"stiff",         2, 0, 's'},
    {"mask",          1, 0, 'M'},
    {"mahnev",        1, 0, 'n'},
    {"mahref",        1, 0, 'x'},
    {"verbose",       1, 0, 'v'},
    {"help",          0, 0, 'h'},
    {"threshold",     1, 0, 't'},
    {0,               0, 0, 0}
  };
  int c, option_index;
  char *strParams = "i:o:Iae:p:r:bf::EA:Plc:T:C:H::s::M:n:R:m:x:vht:";
  CLParameters *clp;

  /* We initialize some variables in order to do parameter recognition */
  option_index = 0;
  clp = (CLParameters *)smalloc (sizeof (*clp));
  clp->infile                = NULL;
  clp->outfile               = NULL;
  clp->maskString            = NULL;
  clp->xfile                 = NULL;
  clp->nEvec                 = 0;
  clp->nProj                 = 0;
  clp->nRms                  = 0;
  clp->nFluc                 = 0;
  clp->nAnim                 = 0;
  clp->nEvecHinge            = 0;
  clp->nEvecCollectivity     = 0;
  clp->operationsAsked       = 0;
  clp->nev                   = 0;
  clp->temperature           = 0.0;
  clp->infoRequested         = FALSE;
  clp->averageRequested      = FALSE;
  clp->evalsRequested        = FALSE;
  clp->evecsRequested        = FALSE;
  clp->projRequested         = FALSE;
  clp->animRequested         = FALSE;
  clp->bfactorsRequested     = FALSE;
  clp->flucRequested         = FALSE;
  clp->lindemannRequested    = FALSE;
  clp->forceCteRequested     = FALSE;
  clp->hingeRequested        = FALSE;
  clp->stiffRequested        = FALSE;
  clp->pdbFormatRequested    = FALSE;
  clp->rmsRequested          = FALSE;
  clp->collectivityRequested = FALSE;
  clp->mahalanobisRequested  = FALSE;
  clp->verbose               = FALSE;
  threshold             = 0.0;

  /* Loop for recognizing input parameters */
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
    case 'I':
      clp->infoRequested = TRUE;
      clp->operationsAsked++;
      break;
    case 'a':
      clp->averageRequested = TRUE;
      clp->operationsAsked++;
      break;
    case 'e':
      clp->evecsRequested = TRUE;
      clp->nEvec = atoi (optarg);
      clp->operationsAsked++;
      break;
    case 'p':
      clp->projRequested = TRUE;
      clp->nProj = atoi (optarg);
      clp->operationsAsked++;
      break;
    case 'r':
      clp->rmsRequested = TRUE;
      clp->nRms = atoi (optarg);
      clp->operationsAsked++;
      break;
    case 'b':
      clp->bfactorsRequested = TRUE;
      break;
    case 'f':
      clp->flucRequested = TRUE;
      if (optarg)
        clp->nFluc = atoi (optarg);
      else
        clp->nFluc = 0;
      clp->operationsAsked++;
      break;
    case 'E':
      clp->evalsRequested = TRUE;
      clp->operationsAsked++;
      break;
    case 'A':
      clp->animRequested = TRUE;
      clp->nAnim = atoi (optarg);
      clp->operationsAsked++;
      break;
    case 'P':
      clp->pdbFormatRequested = TRUE;
      break;
    case 'l':
      clp->lindemannRequested = TRUE;
      clp->operationsAsked++;
      break;
    case 'c':
      clp->forceCteRequested = TRUE;
      clp->temperature = (float)atof (optarg);
      clp->operationsAsked++;
      break;
    case 'T':
      clp->temperature = (float)atof (optarg);
      break;
    case 'C':
      clp->collectivityRequested = TRUE;
      if (optarg)
        clp->nEvecCollectivity = atoi (optarg);
      else
        clp->nEvecCollectivity = 0;
      clp->operationsAsked++;
      break;
    case 'H':
      clp->hingeRequested = TRUE;
      clp->operationsAsked++;
      if (optarg)
        clp->nEvecHinge = atoi (optarg);
      else
        clp->nEvecHinge = 0;
      break;
    case 's':
      clp->stiffRequested = TRUE;
      clp->operationsAsked++;
      if (optarg)
        clp->nEvecStiff = atoi (optarg);
      else
        clp->nEvecStiff = 0;
      break;
    case 'M':
      clp->maskString = copyString (optarg);
      break;
    case 'n':
      clp->nev = atoi (optarg);
      break;
    case 'x':
      clp->mahalanobisRequested = TRUE;
      clp->xfile = copyString (optarg);
      clp->operationsAsked++;
      break;
    case 'v':
      clp->verbose = TRUE;
      break;
    case 'h':
      helptext ();
      exit (0);
      break;
    case 't':
      threshold = atof (optarg);
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
  if (clp->infile != NULL);
    sfree (clp->infile);
  if (clp->outfile != NULL);
    sfree (clp->outfile);
  if (clp->maskString != NULL);
    sfree (clp->maskString);
  if (clp->xfile != NULL);
    sfree (clp->xfile);
  sfree (clp);
}

void checkCLP (CLParameters *clp, PCZFILE *pczf) {
  /* We check that no more than one analysis is requested at one time */
  if (clp->operationsAsked != 1) {
    fprintf (stderr, "You must select one and only one operation to perform\n");
    exit (1);
  }
  
  /* Now we can check that we have correct parameters */
  if (clp->nev == 0) {
    clp->nev = pczf->nvecs;
  }
  if (clp->evecsRequested == TRUE && (clp->nEvec <= 0 || clp->nEvec > pczf->nvecs)) {
	    fprintf (stderr, "Error: There are not so many eigenvectors\n");
	    exit (-1);
  }
  if (clp->projRequested == TRUE && (clp->nProj <= 0 || clp->nProj > pczf->nvecs)) {
	    fprintf (stderr, "Error: There are not so many eigenvectors\n");
	    exit (-1);
  }
  if (clp->flucRequested == TRUE && (clp->nFluc < 0 || clp->nFluc > pczf->nvecs)) {
	    fprintf (stderr, "Error: There are not so many eigenvectors\n");
	    exit (-1);
  }
  if (clp->hingeRequested == TRUE && (clp->nEvecHinge < 0 || clp->nEvecHinge > pczf->nvecs)) {
	    fprintf (stderr, "Error: There are not so many eigenvectors\n");
	    exit (-1);
  }
  if (clp->animRequested == TRUE && (clp->nAnim <= 0 || clp->nAnim > pczf->nvecs)) {
	    fprintf (stderr, "Error: There are not so many eigenvectors\n");
	    exit (-1);
  }
  if (clp->collectivityRequested == TRUE && (clp->nEvecCollectivity < 0 || clp->nEvecCollectivity > pczf->nvecs)) {
	    fprintf (stderr, "Error: There are not so many eigenvectors\n");
	    exit (-1);
  }
  if (clp->mahalanobisRequested == TRUE && (clp->nev <=0 || clp->nev > pczf->nvecs)) {
	    fprintf (stderr, "Error: There are not so many eigenvectors\n");
	    exit (-1);
  }
  if (clp->stiffRequested == TRUE && (clp->temperature <= 0.0)) {
      fprintf (stderr, "Error: Cannot compute stiffness without a temperature\n");
      exit (-1);
  }

  if (clp->nRms < 0 || clp->nRms > pczf->nframes) {
    fprintf (stderr, "Error: There are not so many frames\n");
    exit (-1);
  }
  if (clp->mahalanobisRequested && clp->xfile == NULL) {
    fprintf (stderr, "Error: reference structures file needed\n");
    exit (-1);
  }

  /* We show the user which parameters we readed */
  if (clp->verbose) {
    printf ("Parameters readed:\n");
    printf ("\tInput file: %s\n", clp->infile);
    printf ("\tOutput file: %s\n", clp->outfile);
    if (clp->maskString != FALSE)
      printf ("\tMask: %s\n", clp->maskString);
    if (clp->infoRequested != FALSE)
      printf ("\tInfo: %i\n", clp->infoRequested);
    if (clp->averageRequested != FALSE)
      printf ("\tAverage: %i\n", clp->averageRequested);
    if (clp->nEvec != 0)
      printf ("\tEigenvector: %i\n", clp->nEvec);
    if (clp->nProj != 0)
      printf ("\tProjection: %i\n", clp->nProj);
    if (clp->rmsRequested != 0)
      printf ("\tRMS: %i\n", clp->nRms);
    if (clp->flucRequested != 0)
      printf ("\tFluctuation: %i\n", clp->nFluc);
    if (clp->bfactorsRequested != FALSE)
      printf ("\tB-Factors: %i\n", clp->bfactorsRequested);
    if (clp->hingeRequested != FALSE)
      printf ("\tHinge points: %i (evec %i)\n", clp->hingeRequested, clp->nEvecHinge);
    if (clp->lindemannRequested != FALSE)
      printf ("\tLindemann: %i\n", clp->lindemannRequested);
    if (clp->forceCteRequested != FALSE)
      printf ("\tForce constants: %i\n", clp->forceCteRequested);
    if (clp->collectivityRequested != FALSE)
      printf ("\tCollectivity: %i (evec %i)\n", clp->collectivityRequested, clp->nEvecCollectivity);
    if (clp->evalsRequested != FALSE)
      printf ("\tEigenvalues: %i\n", clp->evalsRequested);
    if (clp->nAnim != 0)
      printf ("\tAnimation: %i\n", clp->nAnim);
    if (clp->pdbFormatRequested != FALSE)
      printf ("\tPDB format: %i\n", clp->pdbFormatRequested);
    if (clp->mahalanobisRequested != FALSE) {
      printf ("\tMahRef: %s\n", clp->xfile);
      printf ("\tMahNev: %i\n", clp->nev);
    }
  }

}
