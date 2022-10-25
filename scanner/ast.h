#ifndef AST_H_
#define AST_H_

enum type {T_UNDEFINED, T_ATOM, T_RESIDUE, T_CHAIN, T_OR, T_AND, T_NOT, T_GROUPING};
typedef enum type TYPE;

/************************************************/
/* Name list to hold the atom and residue names */
/************************************************/
struct nodeNameList {
	char *name;
	struct nodeNameList *next;
};
typedef struct nodeNameList NODENAMELIST;

struct nameList {
	NODENAMELIST *first;
	NODENAMELIST *last;
	TYPE nameType;
};
typedef struct nameList NAMELIST;

NAMELIST *createNameList (void);
void destroyNameList (NAMELIST *nl);
NAMELIST *addNameNameList (NAMELIST *nl, char *i);
NAMELIST *addRangeNameList (NAMELIST *nl, char *i, char *j);
NAMELIST *setTypeNameList (NAMELIST *nl, TYPE nameType);
int containedInNameList (NAMELIST *nl, char *n);
void printNameList (NAMELIST *in);

/*****************************************************/
/* Integer list to hold the atom and residue numbers */
/*****************************************************/
struct nodeIntList {
	int val;
	struct nodeIntList *next;
};
typedef struct nodeIntList NODEINTLIST;

struct intList {
	NODEINTLIST *first;
	NODEINTLIST *last;
	TYPE intType;
};
typedef struct intList INTLIST;

INTLIST *createIntList (void);
void destroyIntList (INTLIST *il);
INTLIST *addIntIntList (INTLIST *il, int i);
INTLIST *addRangeIntList (INTLIST *il, int i, int j);
INTLIST *setTypeIntList (INTLIST *il, TYPE intType);
int containedInIntList (INTLIST *il, int n);
void printIntList (INTLIST *il);

/************************************/
/* Abstract Syntax Tree definitions */
/************************************/
struct nodeAst {
	TYPE nodeType;
	struct nodeAst *left;
	struct nodeAst *right;
	INTLIST *il;
	NAMELIST *nl;
};
typedef struct nodeAst NODEAST;

int acceptedAtom (NODEAST *ast, int atomNumber, char *atomName, int residueNumber, char *residueName, char chainName);
void destroyAst (NODEAST *ast);
NODEAST *joinNodes (NODEAST *left, NODEAST *right, INTLIST *il, NAMELIST *nl, TYPE nodeType);
void printAst (NODEAST *ast);
void printNode (NODEAST *ast, int indentation);
void indent (int indentation);

#endif
