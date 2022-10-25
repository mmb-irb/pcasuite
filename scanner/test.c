#include <stdio.h>
#include "parser.h"

int main (int argc, char **argv) {
  char *mask = ":1-15&~@C,CA,N,O&~@*H*&%A-G";
	int accepted;
	
	printf ("Inici: %s\n", mask);
	processString (mask);
	/* yyparse() method sets the ast global variable */
	yyparse ();
	accepted = acceptedAtom (ast, 2, "CB", 10, "ALA", 'A');
	printf ("Accepted atom? %i\n", accepted);
	destroyAst (ast);
	freeFlexResources  ();
	printf ("Final\n");
	
	return 0;
}
