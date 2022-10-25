#ifndef AUX_PARSER_H_
#define AUX_PARSER_H_

#include "../pcz_io.h"
#include "ast.h"

void yy_scan_string(const char *str);
int yyparse (void);

/* ast value is set by the yyparse() method */
extern NODEAST *ast;

/* Defined in Flex file */
void processString (const char *string);
void freeFlexResources (void);

/* Defined in aux_parser.c
 * nSelected is an output value which gives the number of atoms that match the mask
 * returns a list of the selected atoms Ex.:<1, 3, 15, 16, 17> */
int *selectAtoms (PCZFILE *pczf, char *mask, int *nSelected);
void selectAtomMask (PDBATOM *pdba, int nAtoms, char *mask, int *nSelected, int *selection);

#endif
