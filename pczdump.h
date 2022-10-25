#ifndef PCAZIP_H_
#define PCAZIP_H_

#include <stdio.h>
#include "pcz_io.h"
#include "clusterlist.h" 

struct clparameters_t {
  char *infile;
  char *outfile;
  char *maskString;
  char *xfile;

  int nEvec;
  int nProj;
  int nRms;
  int nFluc;
  int nAnim;
  int nEvecCollectivity;
  int nEvecHinge;
  int nEvecStiff;
  int nev;
  float temperature;

  int infoRequested;
  int averageRequested;
  int rmsRequested;
  int bfactorsRequested;
  int flucRequested;
  int evalsRequested;
  int evecsRequested;
  int projRequested;
  int animRequested;
  int pdbFormatRequested;
  int lindemannRequested;
  int forceCteRequested;
  int collectivityRequested;
  int hingeRequested;
  int stiffRequested;
  int mahalanobisRequested;
  int operationsAsked;
  int verbose;
};
typedef struct clparameters_t CLParameters;

struct dist_t {
	int residue;
	float distance;
};

void dumpevec (PCZFILE *pczf, FILE *fout, int nEvec, int *atomList, int nAtoms);
void dumpevals (PCZFILE *pczf, FILE *fout);
void dumpproj (PCZFILE *pczf, FILE *fout, int ip);
void dumpavgx (PCZFILE *pczf, FILE *fout, int *atomList, int nAtoms);
void dumpfluc (PCZFILE *pczf, FILE *fout, int iv, int bfactorsRequested, int *atomList, int nAtoms);
void dumpflucpdb (PCZFILE *pczf, FILE *fout, int iv, int bfactorsRequested, int *atomList, int nAtoms);
void dumprms (PCZFILE *pczf, FILE *fout, int ir);
void dumpanimx (PCZFILE *pczf, FILE *fout, int iv, int *atomList, int nAtoms);
void dumpavgpdb(PCZFILE *pczf, FILE *fout, int *atomList, int nAtoms/*, char **pdbline, int nr*/);
void dumpanimpdb (PCZFILE *pczf, FILE *fout, int iv, int *atomList, int nAtoms/*, char **pdbline, int nr*/);
void dumpinfo (PCZFILE *pczf, FILE *fout);
void dumpcartfluc (PCZFILE *pczf, FILE *fout, int bfactorsRequested);
void dumplindemann (PCZFILE *pczf, FILE *fout, int *atomList, int nAtoms);
void dumpforcecte (PCZFILE *pczf, FILE *fout, float temperature);
void dumpstiffness (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList, float temperature);
void dumphinge (PCZFILE *pczf, FILE *fout, int nEvec);
void helptext (void);
/*
int **getNeighbourList (PCZFILE *pczf, int *atomList, int nResidues);
int compareDistances (const void *m1, const void *m2);
void getDistances (PCZFILE *pczf, int residue, int *atomList, int nResidues, struct dist_t *distances);
*/
float distance (float *xyz, int i, int j);
float clusterDistance (ClusterList *cl, int clusterA, int clusterB, float **drij);
float hierarchicalClusterDistance (ClusterList *cl, int clusterA, int clusterB, float **drij);
float meanDistanceElemToCluster (ClusterList *cl, int clusterId, float **drij, int elem);
void calchinge_domains (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList);
void getAverageDistances (float *xyz, int nResidues, float *distances);
void dumpcollectivity (PCZFILE *pczf, FILE *fout, int nEvecCollectivity);
void calchinge_distances (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList);
void calchinge_fluc (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList);
void dumpevecpdb(PCZFILE *pczf, FILE *fout, int nEvec, int *atomList, int nAtoms);

/*void calchinge_neighbours (PCZFILE *pczf, FILE *fout, int nEvec, int nResidues, int *atomList);*/
float *meanNeighboursDistanceFluctuations (PCZFILE *pczf, int nEvec, int nResidues, int *atomList, int nNeighbours);
struct dist_t **getNeighboursExtended (float *xyz, int nResidues, int nNeighbours);
struct dist_t *getDistances (float *xyz, int nResidues, int nNeighbours, int residue);
void addDistanceList (struct dist_t *distances, int nNeighbours, int residue, float distance);
float getDistance (float *xyz, int residue1, int residue2);
int **getNeighbours (float *xyz, int nResidues, int nNeighbours);
float getMeanDistance (float *xyz, int residue, int *neighbours, int nNeighbours);

void dumpmahalanobis (PCZFILE *pczf, FILE *fout, int nev, char *xfile, int *atomList, int nAtoms);

CLParameters *readCLParameters (int argc, char **argv);
void freeCLParameters (CLParameters *clp);
void checkCLP (CLParameters *clp, PCZFILE *pczf);

/* roundf declaration. Not found in math.h :-? */
float roundf(float x);

#endif
