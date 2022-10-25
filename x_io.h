#ifndef X_IO_H_
#define X_IO_H_

#include "traj_io.h"

void xopen (TRAJFILE *xf);
int xsnap (TRAJFILE *xf, int is, float *buffer);
int xFormatOK (TRAJFILE *xf);
int checkLineOfNumbers (char *line);
int checkLinePrintable (char *line);
void xclose (TRAJFILE *xf);

#endif
