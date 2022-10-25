#ifndef BINPOS_IO_H_
#define BINPOS_IO_H_

#include "traj_io.h"

void binposopen (TRAJFILE *xf);
int binpossnap (TRAJFILE *xf, int is, float *buffer);
int binposFormatOK (TRAJFILE *xf);
void binposclose (TRAJFILE *xf);

#endif
