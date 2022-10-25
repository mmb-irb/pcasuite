#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "traj_io.h"
#include "x_io.h"
#include "pcazip.h"
#include "utils.h"

/* This method opens an ASCII formatted trajectory file
 * and retrieves the basic information about the file
 * like the number of atoms, frames... */
void xopen (TRAJFILE *xf) {
  char line[82];
  int i, linesPerFrame;

  xf->file      = NULL;
  xf->title     = NULL;
  xf->nframes   = 0;
  xf->nextframe = 0;
  xf->hasbox    = FALSE;
  for (i = 0; i < 6; i++)
    xf->box[i] = 0.0;

  /* Try to open file. Abort if error occurs */
  xf->file = fopen (xf->filename, "r");
  if (xf->file == NULL) {
    perror ("Error: cannot open input file");
    exit (1);
  }

  /* Let's check for the box line */
  line[81] = (char)0;
  for (i = 0; i < 1+((3*xf->natoms)/10); i++) {
    fgets (line, 82, xf->file);
    if (strlen(line) <= 0) {
      fprintf (stderr, "Error reading the file %s\n", xf->filename);
      exit (1);
    }
  }
  if ((3*xf->natoms)%10 != 0) {
    fgets (line, 82, xf->file);
    if (strlen (line) <= 0) {
      fprintf (stderr, "Error reading the file %s\n", xf->filename);
      exit (1);
    }
  }

  /* Now read presumpt box line */
  if (fgets (line, 82, xf->file)) {
    if (strlen (line) <= 0) {
      fprintf (stderr, "Error reading the file %s\n", xf->filename);
      exit (1);
    }
    if (strlen (line) < 49)
      xf->hasbox = TRUE;
  }

  /* Now collect title and find number of frames */
  rewind (xf->file);
  fgets (line, 81, xf->file);
  trim (line);
  xf->title = copyString (line);

  xf->nframes = 0;
  linesPerFrame = (xf->hasbox?1:0) + (3*xf->natoms)/10 + (((3*xf->natoms)%10==0)?0:1);

  while (!feof (xf->file)) {
    for (i = 0; i < linesPerFrame; i++) {
      fgets (line, 82, xf->file);
      if (strlen (line) <= 0)
				break;
    }
    if (!feof (xf->file) && strlen (line) > 0)
      xf->nframes++;
  }

  /* Prepare the file to read first frame */
  rewind (xf->file);
  fgets (line, 82, xf->file);
  xf->nextframe = 1;
}

/* This method searchs for a given frame (is) into the file and returns
 * the coordinates stored in that frame in the supplied buffer */
int xsnap (TRAJFILE *xf, int is, float *buffer) {
  int i, j;
  float *x, *y, *z;

  /* is parameter contains the number of the requested frame */
  /* If we request a frame number higher than the total number of frames
  * we return an error */
  if (is > xf->nframes) {
    return -2;
  }

  /* If the requested frame is efore the actual frame pointer, we reset the
  * file pointer and restart the reading */
  if (is < xf->nextframe) {
    rewind (xf->file);
    fgets (xf->title, 82, xf->file);
    xf->nextframe = 1;
  }

  /* We loop over the frames, reading them and storing the contents into
   * buffer parameter. We finish when we have just readed the requested
   * frame */
  for (i = xf->nextframe; i <= is; i++) {
    x = buffer;
    y = buffer + 1;
    z = buffer + 2;
    for (j = 0; j < xf->natoms; j++) {
      fscanf (xf->file, "%8f%8f%8f", x, y, z);
      x += 3;
      y += 3;
      z += 3;
    }
    /* Check this: I only read 3 box numbers, they could be 6 */
    if (xf->hasbox)
      fscanf (xf->file, "%8f%8f%8f", &xf->box[0], &xf->box[1], &xf->box[2]);

    /* Finally, we adjust the frame pointer to the next frame to be read */
    xf->nextframe++;
  }

  return 0;
}

int xFormatOK (TRAJFILE *xf) {
  int formatOK;
  char line[82];
  int i;

  formatOK = TRUE;

  /* Try to open file. Abort if error occurs */
  xf->file = fopen (xf->filename, "r");
  if (xf->file == NULL) {
    perror ("Error: cannot open input file");
    exit (1);
  }

  /* Let's check for the box line */
  line[81] = (char)0;
  fgets (line, 82, xf->file);
  if (!checkLinePrintable (line)) {
    formatOK = FALSE;
  }
  for (i = 0; i < ((3*xf->natoms)/10); i++) {
    fgets (line, 82, xf->file);
    if (strlen(line) <= 0) {
      fprintf (stderr, "Error reading the file %s\n", xf->filename);
      exit (1);
    }
    if (!checkLineOfNumbers (line)) {
      formatOK = FALSE;
      i = ((3*xf->natoms)/10); /* This is the exit condition */
    }
  }
  if ((3*xf->natoms)%10 != 0) {
    fgets (line, 82, xf->file);
    if (strlen (line) <= 0) {
      fprintf (stderr, "Error reading the file %s\n", xf->filename);
      exit (1);
    }
    if (!checkLineOfNumbers (line)) {
      formatOK = FALSE;
    }
  }

  /* Now read presumpt box line */
  if (fgets (line, 82, xf->file)) {
    if (strlen (line) <= 0) {
      fprintf (stderr, "Error reading the file %s\n", xf->filename);
      exit (1);
    }
    if (!checkLineOfNumbers (line)) {
      formatOK = FALSE;
    }
  }

  fclose (xf->file);

  return formatOK;
}

void xclose (TRAJFILE *xf) {
  fclose (xf->file);
}

/* This method checks if all the characters in the line are numbers */
int checkLineOfNumbers (char *line) {
  char *ptr;
  int isOK;
  
  isOK = TRUE;
  ptr = line;
  while (*ptr != (char)0 && isOK) {
  	isOK = isOK && (isdigit (*ptr) || isspace (*ptr) || *ptr == '.' || *ptr == '-');
  	ptr++;
  }
  
  return isOK;
}

/* This method checks if all the characters in the line are printble */
int checkLinePrintable (char *line) {
  char *ptr;
  int isOK;
  
  isOK = TRUE;
  ptr = line;
  while (*ptr != (char)0 && isOK) {
  	isOK = isOK && (isprint (*ptr) || isspace (*ptr));
  	ptr++;
  }
  
  return isOK;
}
