#ifndef CONSTANTS_H_
#define CONSTANTS_H_

/* Constants */
/*************/
#define FALSE 0
#define TRUE !FALSE

#define MAX_LINE_LENGTH 82

#define NONE_RMS    -1
#define STD_RMS      0
#define GAUSSIAN_RMS 1

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define CODE_VERSION "2.1.0"

#define VERSION "PCZ4"

/* Value under we consider data as zero */
#define EPSILON 0.0000001

/* r -> J/K.mol */
#define R_CTE 8.314

/* jules to Kcal */
#define F_CTE (4.184 * 1000)

/* Structures */
/**************/
struct pdbatom_t {
	int  atomNumber;
	int  residueNumber;
	char residueName[4];
	char atomName[5];
  char chain;
};
typedef struct pdbatom_t PDBATOM;

/* This is the contant used through the programs to store the verbosity level
 * desired by the user */
extern int verbose;

#endif
