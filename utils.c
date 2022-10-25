#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "utils.h"

/* File existence test. It doesn't check if it can be readed */
int fexist (char *filename) {
  struct stat buffer;
  
  if (stat (filename, &buffer) == 0)
    return 1; /* Stat successful, file exists */
  else
    return 0; /* Stat unsuccessful, file cannot be read */
}

/* This method copies the input string in a newly allocated memory junk */
char *copyString (char *str) {
  char *result;
  size_t length;

  length = strlen (str);
  result = (char *) smalloc ((length + 1) * sizeof (char));
  strcpy (result, str);

  return result;
}

/* This method returns the minimum value of the vector elements */
float minval (float *v, int nElems) {
  float val;
  int i;

  val = v[0];
  for (i = 1; i < nElems; i++) {
    val = (val > v[i])? v[i] : val;
  }

  return val;
}

/* This method returns the maximum value of the vector elements */
float maxval (float *v, int nElems) {
  float val;
  int i;

  val = v[0];
  for (i = 1; i < nElems; i++) {
    val = (val < v[i])? v[i] : val;
  }

  return val;
}

/* This method strips the blank characters at the beginning and the end of a
 * string and leaves the result at the beginning of the original string */
char *trim (char *b) {
  char *e, *aux;

  aux = b;
  e = strrchr(b, '\0'); /* Find the final null */
  while (aux < e && isspace (*aux)) /* Scan forward */
    ++aux;
  while (e > aux && isspace (*(e-1))) /* scan back from end */
    --e;
  *e = '\0'; /* terminate new string */

  memmove (b, aux, strlen (aux) + 1);

  return b;
}

/* Output of floats in (10f8.3) format */
void writeVector (float *vector, int nelems, FILE *fout) {
  int colElems, i;

  colElems = 0;
  for (i = 0; i < nelems; i++) {
    if (colElems == 10) {
      colElems = 0;
      fprintf (fout, "\n");
    }
    fprintf (fout, "%8.3f", vector[i]);
    colElems++;
  }
  fprintf (fout, "\n");
}

/* Output of doubles in (10f8.3) format */
void writeVectorD (double *vector, int nelems, FILE *fout) {
  int colElems, i;

  colElems = 0;
  for (i = 0; i < nelems; i++) {
    if (colElems == 10) {
      colElems = 0;
      fprintf (fout, "\n");
    }
    fprintf (fout, "%8.3f", vector[i]);
    colElems++;
  }
  fprintf (fout, "\n");
}

/* Output of floats in (f10.6) format */
void writeColumnVector (float *vector, int nelems, FILE *fout) {
  int i;

  for (i = 0; i < nelems; i++) {
    fprintf (fout, "%10.6f\n", vector[i]);
  }
}

/* Output of floats in CSV format */
void writeCSV (float *vector, int nelems, FILE *fout) {
  int i;

  fprintf (fout, "%8.3f", vector[0]);
  for (i = 1; i < nelems; i++) {
    fprintf (fout, ",%8.3f", vector[i]);
  }
  fprintf (fout, "\n");
}

/* Output of floats in CSV format and Scientific notation */
void writeCSVSci (float *vector, int nelems, FILE *fout) {
  int i;

  fprintf (fout, "%.3e", vector[0]);
  for (i = 1; i < nelems; i++) {
    fprintf (fout, ",%.3e", vector[i]);
  }
  fprintf (fout, "\n");
}

/* This method generates and returns a new vector with the coordinates
 * for the requested atoms copied. This method suppose that the input data
 * are atom coordinates and copies the *atoms*, this is the three coordinates
 * for each atom.
 * xyz is the source vector,
 * atomList is a boolean vector containing 1 for elements that should be
 *     copied and 0 for elements than must not be copied
 * nAtoms is the number of elements in each list */
float *packAtoms (float *xyz, int *atomList, int nAtoms) {
  int i;
  float *selection;

  selection = (float *)smalloc (nAtoms * 3 * sizeof (float));
#pragma omp parallel for schedule(static) shared(selection, xyz, atomList) private(i) firstprivate(nAtoms)
  for (i = 0; i < nAtoms; i++) {
    selection[i*3  ] = xyz[atomList[i]*3  ];
    selection[i*3+1] = xyz[atomList[i]*3+1];
    selection[i*3+2] = xyz[atomList[i]*3+2];
  }

  return selection;
}

/* This method generates and returns a new vector with the coordinates
 * for the requested values copied.
 * xyz is the source vector,
 * atomList is a boolean vector containing 1 for elements that should be
 *     copied and 0 for elements than must not be copied
 * nAtoms is the number of elements in each list */
float *packValues (float *values, int *atomList, int nAtoms) {
  int i;
  float *selection;

  selection = (float *)smalloc (nAtoms * sizeof (float));
#pragma omp parallel for schedule(static) shared(selection, values, atomList) private(i) firstprivate(nAtoms)
  for (i = 0; i < nAtoms; i++) {
    selection[i] = values[atomList[i]];
  }

  return selection;
}

/* This method computes the mean value and the variance of a list of integers.
 * The results are stored in the mean and variance out parameters */
void computeMeanVariance (int *values, int nValues, float *mean, float *variance) {
  int i;

  *mean = 0.0;
  *variance = 0.0;
  for (i = 0; i < nValues; i++) {
    *mean += (float)values[i];
    *variance += (float)values[i] * (float)values[i];
  }
  *mean /= (float)nValues;
  *variance = (*variance / (float)nValues) - (*mean * *mean);
}

/* Quick test to check for the endianness of the executing machine */
int am_big_endian (void) {
   long one = 1;
   return !(*((char *)(&one)));
}

/* Method to change the endianness of a float value between little endian
 * and big endian */
float swap_endiannessf (float f) {
  float b;
  unsigned int ui, ui2;

  ui = *(unsigned int *)(void *)&f;
  ui2 = swap_endianness (ui);

  b = *(float *)(void *)&ui2;

  return b;
}

/* Method to change the endianness of an integer value between little endian
 * and big endian */
unsigned int swap_endianness (unsigned int a) {
  unsigned int a1, a2, a3, a4, b;

  a1 = a & 0xFF000000;
  a2 = a & 0x00FF0000;
  a3 = a & 0x0000FF00;
  a4 = a & 0x000000FF;

  a1 = (a1 >> (8+8+8)) & 0x000000FF;
  a2 = (a2 >> 8) & 0x0000FF00;
  a3 = a3 << 8;
  a4 = a4 << (8+8+8);

  b = a1 | a2 | a3 | a4;

  return b;
}

/* Rotates a string one character to the left, having the discarded character
 * enter into the right side "ABCD"-> "BCDA" */
char *rotateLeft(char *name, int length) {
	int i;
	char carry;

	carry = name[0];
	for (i = 0; i < length-1; i++) {
		name[i] = name[i+1];
	}
	name[length-1] = carry;

	return name;
}

/* Rotates a string one character to the right, having the discarded character
 * enter into the left side "ABCD"-> "DABC" */
char *rotateRight(char *name, int length) {
	int i;
	char carry;

	carry = name[length-1];
	for (i = length-1; i > 0; i--) {
		name[i] = name[i-1];
	}
	name[0] = carry;

	return name;
}

/* Memory allocation routine that performs the basic sanity checks each time
 * memory is requested */
void *smalloc (size_t size) {
	void *ptr;

	ptr = malloc (size);

	if (ptr == NULL) {
		fprintf (stderr, "Not enough memory allocating %li bytes\n", (unsigned long int)size);
		exit (1);
	}

	return ptr;
}

/* Memory deallocation routine that performs the basic sanity checks each time
 * memory is freed */
void sfree (void *ptr) {
	if (ptr != NULL)
	  free (ptr);
}

/* This method outputs a PDB-like file containing the coordinates and
 * b-factors given in the xyz and bfactors parameters, using the names
 * found in the pczf structure */
void dumppdb (PCZFILE *pczf, FILE *fout, float *xyz, float *bfactors, int *atomList, int nAtoms) {
  int i, atomNumber;
  char atomName[5];

  for (i = 0; i < nAtoms; i++) {
    if (atomList != NULL)
      atomNumber = atomList[i];
		else
		  atomNumber = i;
    if (strlen(pczf->atomNames[atomNumber].atomName) == 4) {
      sprintf (atomName, "%4s", pczf->atomNames[atomNumber].atomName);
      rotateRight(atomName, strlen(atomName));
    } else {
      sprintf (atomName, " %-3s", pczf->atomNames[atomNumber].atomName);
    }
    if (bfactors != NULL) {
      fprintf (fout, "ATOM  %5i %4.4s %3.3s %c%4i    %8.3f%8.3f%8.3f      %6.2f\n",
        pczf->atomNames[atomNumber].atomNumber,
        atomName,
        pczf->atomNames[atomNumber].residueName,
        pczf->atomNames[atomNumber].chain,
        pczf->atomNames[atomNumber].residueNumber,
        xyz[3*i],
        xyz[3*i+1],
        xyz[3*i+2],
        bfactors[i]
        );
     } else {
      fprintf (fout, "ATOM  %5i %4.4s %3.3s %c%4i    %8.3f%8.3f%8.3f\n",
        pczf->atomNames[atomNumber].atomNumber,
        atomName,
        pczf->atomNames[atomNumber].residueName,
        pczf->atomNames[atomNumber].chain,
        pczf->atomNames[atomNumber].residueNumber,
        xyz[3*i],
        xyz[3*i+1],
        xyz[3*i+2]
        );
     }
  }
}

float *matmul (float **matrix, float *vector, int rows, int columns) {
  float *result;
  int i, j;
  
  result = (float *)smalloc (rows * sizeof (float));
  memset (result, 0, rows * sizeof (float));
  
#pragma omp parallel for schedule(static) shared (result, matrix, vector, rows, columns) private (i, j)
  for (i = 0; i < rows; i++) {
    for (j = 0; j < columns; j++) {
      result[i] += matrix[i][j] * vector[j];
    }
  }
  
  return result;
}

long int getTime (void) {
  struct timeval  tp; 
  long int result;
  
  gettimeofday(&tp, NULL); 
  result = tp.tv_sec * 1000000 + tp.tv_usec;
  
  return result;
}

/* Method to read a PDB file */
void pdbread (char *pdbfile, int natom, PDBATOM *pdba) {
  FILE *pdbf;
  char line[82];
  char str[6];
  int atomNumber	;

  atomNumber = 0;
  pdbf = fopen (pdbfile, "r");
  fgets (line, 82, pdbf);
  while (!feof (pdbf) && (natom > atomNumber)) {
    if (strncmp (line, "ATOM  ", 6) == 0) {
      strncpy (str, line+6, 5);
      str[5] = (char)0;
      pdba[atomNumber].atomNumber = atoi (str);

      strncpy (str, line+22, 4);
      str[4] = (char)0;
      pdba[atomNumber].residueNumber = atoi (str);

      strncpy (pdba[atomNumber].atomName, line+12, 4);
      pdba[atomNumber].atomName[4] = (char)0;

      strncpy (pdba[atomNumber].residueName, line+17, 3);
      pdba[atomNumber].residueName[3] = (char)0;

      pdba[atomNumber].chain = line[21];
      
      atomNumber++;
    }
    fgets (line, 82, pdbf);
  }
  fclose (pdbf);
}

int addToVector (int *vector, int elem, int next) {
  int idx, *found;
  
  /* First we search for the existence of the element */
  found = bsearch (&elem, vector, next, sizeof (*vector), cmpintegers);
  
  if (found != NULL) {
    /* We found the element, no changes needed to the vector */
    idx = next;
  } else {
    /* The element is not there, so we add the element at the end */
    vector[next] = elem;
    idx = next + 1;
  }
  
  /* Finally we order the vector to speedup the access */
  qsort (vector, idx, sizeof (*vector), cmpintegers); 
  
  return idx;
}

int cmpintegers(const void *p1, const void *p2) {
  return *(int *)p1 - *(int *)p2;
}
