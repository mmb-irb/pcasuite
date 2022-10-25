#ifndef TRAJ_IO_H_
#define TRAJ_IO_H_

#include <stdio.h>

typedef struct trajfile_t TRAJFILE;
typedef void (*trjopen_t)(TRAJFILE *);
typedef int  (*trjsnap_t)(TRAJFILE *, int, float *);
typedef int  (*trjFormatOK_t)(TRAJFILE *);
typedef void (*trjclose_t)(TRAJFILE *);

struct trajfile_t {
  char *filename;
  char *title;
  FILE *file;
  int fd;
  int natoms;
  int nframes;
  int nextframe;
  int hasbox;
  float box[6];
  trjopen_t trjopen;
  trjsnap_t trjsnap;
  trjclose_t trjclose;
};

void trajopen (char *fname, int natoms, TRAJFILE *xf);
void trajclose (TRAJFILE *xf);
void trajfree (TRAJFILE *xf);
int trajsnap (TRAJFILE *xf, int is, float *buffer);
int trajbox (TRAJFILE *xf);
int trajframes (TRAJFILE *xf);
float **trajmemread (TRAJFILE *xf);
int prova (int t);
#endif
