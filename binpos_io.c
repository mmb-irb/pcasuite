#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "pcazip.h"
#include "utils.h"
#include "traj_io.h"

/* This method opens a binpos formatted trajectory file
 * and retrieves the basic information about the file
 * like the number of atoms, frames... */
void binposopen (TRAJFILE *xf) {
  int i;
  size_t peof;

  /* The title will be freed later, so we need to supply memory space */
  xf->title = (char *)malloc (sizeof (char));
  xf->title[0] = (char)0;
  xf->hasbox = FALSE;
  for (i = 0; i < 6; i++)
    xf->box[i] = 0.0;

  /* Try to open file. Abort if error occurs */
  xf->file = fopen (xf->filename, "rb");
  if (xf->file == NULL) {
    perror ("Error: cannot open input file");
    exit (1);
  }

  /* Now get number of atoms and find number of frames */
  fseek (xf->file, 4, SEEK_SET);
  fread (&xf->natoms, 4, 1, xf->file);
  fseek (xf->file, 0, SEEK_END);
  peof = ftell(xf->file);
  xf->nframes = (peof-4) / (((xf->natoms*3)+1) * 4);
  
  /* Prepare the file to read first frame */
  fseek (xf->file, 4, SEEK_SET);
  xf->nextframe = 1;
}

/* This method searchs for a given frame (is) into the file and returns
 * the coordinates stored in that frame in the supplied buffer */
int binpossnap (TRAJFILE *xf, int is, float *buffer) {
  /* is parameter contains the number of the requested frame */
  /* If we request a frame number higher than the total number of frames
   * we return an error */
  if (is > xf->nframes) {
    return -2;
  }

  /* We reset the input pointer to the required register */
  fseek (xf->file, 8+((is-1)*(((xf->natoms*3)+1)*4)), SEEK_SET);

  /* We read the frame and store the contents into buffer parameter. */
  fread (buffer, 4, xf->natoms * 3, xf->file);
  
  /* Finally, we adjust the frame pointer to the next frame to be read */
  xf->nextframe = is + 1;
  
  return 0;
}

/* This method return a boolean value that indicates if the file is in the
 * format recognised by this module */
int binposFormatOK (TRAJFILE *xf) {
  char header[4];
  int formatOK;
  
  formatOK = TRUE;
  xf->file = fopen (xf->filename, "rb");
  if (xf->file == NULL) {
    perror ("Error: cannot open input file");
    exit (1);
  }

  /* Check if the format is the expected */
  fread (header, 1, 4, xf->file);
  if (header[0] != 'f' || header[1] != 'x' || header[2] != 'y' || header[3] != 'z') {
    formatOK = FALSE;
  }
  
  fclose(xf->file);
  
  return formatOK;
}

void binposclose (TRAJFILE *xf) {
  fclose (xf->file);
}
