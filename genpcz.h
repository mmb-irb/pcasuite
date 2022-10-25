#ifndef _GENPCZ_H_
#define _GENPCZ_H_

struct clparameters_t {
  char *evalsFile;
  char *evecsFile;
  char *projFile;
  char *pdbFile;
  char *pczFile;
  char *title;
  int dimensionality;
  int rmsType;
  int verbose;
};
typedef struct clparameters_t CLParameters;

float *readEvals (char *evalsFile, int *nEvals);
float **readEvecs (char *evecsFile, int *nEvecs, float **mean);
float **readProj (char *projFile, int *nEvecs, int *nFrames);
CLParameters *readCLParameters (int argc, char **argv);
void freeCLParameters (CLParameters *clp);
void helptext (void);

#endif
