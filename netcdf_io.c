#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include "pcazip.h"
#include "utils.h"
#include "traj_io.h"

/* This method opens a binpos formatted trajectory file
 * and retrieves the basic information about the file
 * like the number of atoms, frames... */
void netcdfopen (TRAJFILE *xf) {
  int stat, atomId, frameId, coordId;
  size_t titleLength, natoms, nframes;
  int i;

  /* Won't read the box info */
  xf->hasbox = FALSE;
  for (i = 0; i < 6; i++)
    xf->box[i] = 0.0;

  /* Try to open file. Abort if error occurs */
  stat = nc_open (xf->filename, NC_NOWRITE, &xf->fd);
  
  /* Get the variable IDs */
  stat = nc_inq_attlen (xf->fd, NC_GLOBAL, "title", &titleLength);
  stat = nc_inq_dimid (xf->fd, "atom", &atomId);
  stat = nc_inq_dimid (xf->fd, "frame", &frameId);
  stat = nc_inq_varid (xf->fd, "coordinates", &coordId);
  
  /* Get the actual data */
  xf->title = (char *)malloc ((titleLength + 1) * sizeof (char));
  stat = nc_get_att_text (xf->fd, NC_GLOBAL, "title", xf->title);
  xf->title[titleLength] = (char)0;
  stat = nc_inq_dimlen (xf->fd, atomId, &natoms);
  stat = nc_inq_dimlen (xf->fd, frameId, &nframes);
  xf->natoms = (int)natoms;
  xf->nframes = (int)nframes;
  
  /* Prepare for read */
  xf->nextframe = 1;
}

/* This method searchs for a given frame (is) into the file and returns
 * the coordinates stored in that frame in the supplied buffer */
int netcdfsnap (TRAJFILE *xf, int is, float *buffer) {
  int stat, coordId;
  size_t start[3], count[3];
  
  /* is parameter contains the number of the requested frame */
  /* If we request a frame number higher than the total number of frames
   * we return an error */
  if (is > xf->nframes) {
    return -2;
  }

  /* Obtain the variable ID */
  stat = nc_inq_varid (xf->fd, "coordinates", &coordId);
  
  /* Read the data in one chunk */
  start[0] = is - 1; /* Record to read */
  start[1] = 0; /* Start from the beginning of record */
  start[2] = 0; /* Start from the beginning of record */
  count[0] = 1; /* Read only one record */
  count[1] = xf->natoms; /* Read all the record */
  count[2] = 3; /* Read all the record */
  stat = nc_get_vara_float (xf->fd, coordId, start, count, buffer);
  if (stat != NC_NOERR) {
    fprintf (stderr, "Error reading NetCDF file (frame %i)\n", is);
    exit (1);
  }

  /* Finally, we adjust the frame pointer to the next frame to be read */
  xf->nextframe = is + 1;
  
  return 0;
}

/* This method return a boolean value that indicates if the file is in the
 * format recognised by this module */
int netcdfFormatOK (TRAJFILE *xf) {
  int ncid, stat, atomId, frameId, coordId;
  size_t titleLength;
  int formatOK;
    
  formatOK = TRUE;
  
  /* Open the file */
  stat = nc_open (xf->filename, NC_NOWRITE, &ncid);
  if (stat != NC_NOERR) {
    formatOK = FALSE;
  } else {
    /* Check if the format is the expected */
    stat = nc_inq_attlen (ncid, NC_GLOBAL, "title", &titleLength);
    if (stat != NC_NOERR) {
      printf ("Title cannot be read\n");
      formatOK = FALSE;
    }
    stat = nc_inq_dimid (ncid, "atom", &atomId);
    if (stat != NC_NOERR) {
      printf ("Atom cannot be read\n");
      formatOK = FALSE;
    }
    stat = nc_inq_dimid (ncid, "frame", &frameId);
    if (stat != NC_NOERR) {
      printf ("Frame cannot be read\n");
      formatOK = FALSE;
    }
    stat = nc_inq_varid (ncid, "coordinates", &coordId);
    if (stat != NC_NOERR) {
      printf ("Coordinates cannot be read\n");
      formatOK = FALSE;
    }
    
    /* Close the file */
    stat = nc_close (ncid);
  }
  
  return formatOK;
}

void netcdfclose (TRAJFILE *xf) {
  nc_close (xf->fd);
}
