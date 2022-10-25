#include <stdlib.h>
#include <string.h>
#include "aux_parser.h"
#include "../pcz_io.h"
#include "../utils.h"

int *selectAtoms (PCZFILE *pczf, char *mask, int *nSelected) {
	int *selection, i;

	*nSelected = 0;
	selection = (int *)malloc (pczf->natoms * sizeof (int));
	if (pczf->haveAtomNames == FALSE) {
		/* We return a full selection */
		for (i = 0; i < pczf->natoms; i++)
			selection[i] = i;
	} else {
		/* We compute the selection */
		processString (mask);
		yyparse ();
		*nSelected = 0;
		for (i = 0; i < pczf->natoms; i++) {
			if (acceptedAtom (ast,
				pczf->atomNames[i].atomNumber,
				pczf->atomNames[i].atomName,
				pczf->atomNames[i].residueNumber,
				pczf->atomNames[i].residueName,
        pczf->atomNames[i].chain)) {
				selection[*nSelected] = i;
				(*nSelected)++;
			}
		}
		destroyAst (ast);
		freeFlexResources ();
	}

	return selection;
}

void selectAtomMask (PDBATOM *pdba, int nAtoms, char *mask, int *nSelected, int *selection) {
	char atomName[5];
	int i;

	*nSelected = 0;
	/* We compute the selection */
	processString (mask);
	yyparse ();
	*nSelected = 0;
	for (i = 0; i < nAtoms; i++) {
		strcpy (atomName, pdba[i].atomName);
		rotateLeft (atomName, 4);
		trim (atomName);
		if (acceptedAtom (ast,
			pdba[i].atomNumber,
			atomName,
			pdba[i].residueNumber,
			pdba[i].residueName,
      pdba[i].chain)) {
			selection[i*3]   = 1;
			selection[i*3+1] = 1;
			selection[i*3+2] = 1;
			(*nSelected)++;
		} else {
			selection[i*3]   = 0;
			selection[i*3+1] = 0;
			selection[i*3+2] = 0;
		}
	}
	destroyAst (ast);
	freeFlexResources ();
}
