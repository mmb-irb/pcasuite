#ifndef TRAJUTILS_H_
#define TRAJUTILS_H_

void trajfit (float *x, float *ref, int natoms, int nframes);
void trajavg (float *x, float *avg, int natoms, int nframes);

/*#ifdef __xlc__*/
  /* This definition is needed to accomodate the standard way of defining
   * the names in the different compilers and platforms. Because we rely
   * on external libraries, we cannot use the compiler options to control
   * this behaviour, because the external libraries are compiled with
   * default options */
/*  #define matfit_ matfit*/
/*#endif*/

void matfit (int *natoms, float *ref, float *x, float *r,
	     float *v, float *rmsd, int *findmove);

#endif
