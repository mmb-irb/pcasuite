#ifndef NETCDF_IO_H_
#define NETCDF_IO_H_

#include "traj_io.h"

void netcdfopen (TRAJFILE *xf);
int netcdfsnap (TRAJFILE *xf, int is, float *buffer);
int netcdfFormatOK (TRAJFILE *xf);
void netcdfclose (TRAJFILE *xf);

#endif /*NETCDF_IO_H_*/
