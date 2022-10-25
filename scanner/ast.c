#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ast.h"

/**************************************/
/* Abstract Syntax Tree query methods */
/**************************************/
int acceptedAtom (NODEAST *ast, int atomNumber, char *atomName, int residueNumber, char *residueName, char chainName) {
	int accepted;
  char chainId[2];
  
  chainId[0] = chainName;
  chainId[1] = '\0';
	
	switch (ast->nodeType) {
	case T_AND:
	  accepted = acceptedAtom (ast->left, atomNumber, atomName, residueNumber, residueName, chainName) &&
	    acceptedAtom (ast->right, atomNumber, atomName, residueNumber, residueName, chainName);
	  break;
	case T_OR:
	  accepted = acceptedAtom (ast->left, atomNumber, atomName, residueNumber, residueName, chainName) ||
	    acceptedAtom (ast->right, atomNumber, atomName, residueNumber, residueName, chainName);
	  break;
	case T_NOT:
	  accepted = !acceptedAtom (ast->left, atomNumber, atomName, residueNumber, residueName, chainName);
	  break;
	case T_GROUPING:
	  accepted = 0;
	  if (ast->il != NULL) {
			if (ast->il->intType == T_ATOM) {
				accepted = containedInIntList (ast->il, atomNumber);
      } else if (ast->il->intType == T_RESIDUE) {
        accepted = containedInIntList (ast->il, residueNumber);
			}
		} else if (ast->nl != NULL) {
		  if (ast->nl->nameType == T_ATOM) {
				accepted = containedInNameList (ast->nl, atomName);
      } else if (ast->nl->nameType == T_RESIDUE) {
        accepted = containedInNameList (ast->nl, residueName);
      } else if (ast->nl->nameType == T_CHAIN) {
        accepted = containedInNameList (ast->nl, chainId);
			}
		}
	  break;
	default:
	  fprintf (stderr, "Error: Wrong AST node type\n");
	  exit(3);
	  break;
	}

	return accepted;
}

/********************************/
/* Abstract Syntax Tree methods */
/********************************/
NODEAST *joinNodes (NODEAST *left, NODEAST *right, INTLIST *il, NAMELIST *nl, TYPE nodeType) {
	NODEAST *node;

	node = (NODEAST *)malloc (sizeof (NODEAST));
	node->nodeType = nodeType;
	node->left = left;
	node->right = right;
	node->il = il;
	node->nl = nl;

	return node;
}

void destroyAst (NODEAST *ast) {
	if (ast == NULL)
	  return;
	  
	if (ast->left != NULL)
	  destroyAst (ast->left);
	if (ast->right != NULL)
	  destroyAst (ast->right);
	if (ast->il != NULL)
	  destroyIntList (ast->il);
	if (ast->nl != NULL)
	  destroyNameList (ast->nl);
	free (ast);
}

void printAst (NODEAST *ast) {
	if (ast == NULL)
	  return;

	printNode (ast, 0);
}

void printNode (NODEAST *ast, int indentation) {
	indent (indentation);
	switch (ast->nodeType) {
	case T_OR:
	  printf ("OR: ");
	  break;
	case T_AND:
	  printf ("AND: ");
	  break;
	case T_NOT:
	  printf ("NOT: ");
	  break;
	case T_GROUPING:
	  printf ("GROUPING: ");
	  break;
	default:
	  printf ("No Type: ");
	  break;
	}
	if (ast->il != NULL)
	  printIntList (ast->il);
	if (ast->nl != NULL)
	  printNameList (ast->nl);
	printf ("\n");
	
	if (ast->left != NULL)
		printNode (ast->left, indentation + 1);
	if (ast->right != NULL)
		printNode (ast->right, indentation + 1);
}

void indent (int indentation) {
	int i;
	
	for (i = 0; i < indentation * 2; i++)
	  printf (" ");
}

/************************************************/
/* Name list to hold the atom and residue names */
/************************************************/
NAMELIST *createNameList (void) {
	NAMELIST *nl;

	nl = (NAMELIST *)malloc (sizeof (NAMELIST));
	nl->nameType = T_UNDEFINED;
	nl->first = NULL;
	nl->last = NULL;

	return nl;
}

void destroyNameList (NAMELIST *nl) {
	NODENAMELIST *node;

	if (nl == NULL)
	  return;

	node = nl->first;
	while (node != NULL) {
		node = node->next;
		free (nl->first->name);
		free (nl->first);
		nl->first = node;
	}
	free (nl);
}

NAMELIST *addNameNameList (NAMELIST *nl, char *i) {
	NODENAMELIST *node;

	if (nl == NULL)
	  nl = createNameList();

	node = (NODENAMELIST *)malloc (sizeof (NODENAMELIST));
	node->name = (char *)malloc ((strlen (i) + 1) * sizeof (char));;
	strcpy (node->name, i);
	node->next = NULL;

	if (nl->first == NULL) {
		nl->first = node;
	} else {
	  nl->last->next = node;
	}
	nl->last = node;

	return nl;
}

NAMELIST *addRangeNameList (NAMELIST *nl, char *i, char *j) {
  char r[2];
  r[1] = '\0';
  
  if (i[0] > j[0]) {
    fprintf (stderr, "Error in range: %c-%c\n", i[0], j[0]);
    exit(3);
  }
  
  for (r[0] = i[0]; r[0] <= j[0]; r[0]++) {
    nl = addNameNameList (nl, r);
  }

  return nl;

}

NAMELIST *setTypeNameList (NAMELIST *nl, TYPE nameType) {
	nl->nameType = nameType;
	
	return nl;
}

int containedInNameList (NAMELIST *nl, char *n) {
	NODENAMELIST *node;
	int found, length;

	found = 0;
	node = nl->first;
	while (node != NULL && !found) {
		length = strlen (node->name);
		if (node->name[0] == '*' && node->name[length - 1] != '*') {
		  if (length - 1 > strlen (n)) {
				found = 0;
			} else {
				found = strcmp (&node->name[1], n + strlen (n) - length + 1) == 0;
			}
		} else if (node->name[0] == '*' && node->name[length - 1] == '*') {
			node->name[length - 1] = (char)0;
			found = strstr (n, &node->name[1]) != NULL;
			node->name[length - 1] = '*';
		} else if (node->name[0] != '*' && node->name[length - 1] == '*') {
			found = strncmp (node->name, n, length - 1) == 0;
		} else {
			found = strcmp (node->name, n) == 0;
		}
		node = node->next;
	}
	
	return found;
}

void printNameList (NAMELIST *nl) {
	NODENAMELIST *node;

	if (nl == NULL)
		return;

	if (nl->nameType == T_ATOM) {
		printf ("Atoms: ");
  } else if (nl->nameType == T_RESIDUE) {
    printf ("Residues: ");
  } else if (nl->nameType == T_CHAIN) {
    printf ("Chains: ");
	} else {
		printf ("Generic names: ");
	}
	node = nl->first;
	while (node != NULL) {
		printf ("%s, ", node->name);
		node = node->next;
	}
}

/*****************************************************/
/* Integer list to hold the atom and residue numbers */
/*****************************************************/
INTLIST *createIntList (void) {
	INTLIST *il;
	
	il = (INTLIST *)malloc (sizeof (INTLIST));
	il->intType = T_UNDEFINED;
	il->first = NULL;
	il->last = NULL;
	
	return il;
}

void destroyIntList (INTLIST *il) {
	NODEINTLIST *node;
	
	if (il == NULL)
	  return;
	  
	node = il->first;
	while (node != NULL) {
		node = node->next;
		free (il->first);
		il->first = node;
	}
	free (il);
}

INTLIST *addIntIntList (INTLIST *il, int i) {
	NODEINTLIST *node;
	
	if (il == NULL)
	  il = createIntList();
	
	node = (NODEINTLIST *)malloc (sizeof (NODEINTLIST));
	node->val = i;
	node->next = NULL;

	if (il->first == NULL) {
		il->first = node;
	} else {
	  il->last->next = node;
	}
	il->last = node;

	return il;
}

INTLIST *addRangeIntList (INTLIST *il, int i, int j) {
  int r;

  if (i > j) {
    fprintf (stderr, "Error in range: %i-%i\n", i, j);
    exit(3);
  }
  
  for (r = i; r <= j; r++) {
    il = addIntIntList (il, r);
  }

  return il;

}

INTLIST *setTypeIntList (INTLIST *il, TYPE intType) {
	il->intType = intType;
	
	return il;
}

int containedInIntList (INTLIST *il, int n) {
	NODEINTLIST *node;
	int found;
	
	found = 0;
	node = il->first;
	while (node != NULL && !found) {
		found = node->val == n;
		node = node->next;
	}
	
	return found;
}

void printIntList (INTLIST *il) {
	NODEINTLIST *node;
	
	if (il == NULL)
		return;

	if (il->intType == T_ATOM) {
		printf ("Atoms: ");
	} else if (il->intType == T_RESIDUE) {
		printf ("Residues: ");
	} else {
		printf ("Generic ints: ");
	}
	node = il->first;
	while (node != NULL) {
		printf ("%i, ", node->val);
		node = node->next;
	}
}
