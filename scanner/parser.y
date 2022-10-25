%{
  #include <stdio.h>
  #include "ast.h"
  int yylex (void);
  void yyerror (char const *);
  NODEAST *ast;
%}

%union {
	int num;
	char *name;
	char *chainid;
	INTLIST *il;
	NAMELIST *nl;
	NODEAST *ast;
};

%type <il> numlist
%type <ast> atomspec residuespec chainspec unitspec expression root
%type <nl> namelist chainlist

%token <num> NUM
%token <name> NAME
%token RESIDUE
%token ATOM
%token CHAIN
%token OPEN
%token CLOSE
%token COMMA
%left AND OR
%left NOT
%right RANGE

%%

root: expression {$$ = $1; ast = $$; /*printAst ($$);*/}
	;

expression: unitspec {$$ = $1;}
	| NOT expression {$$ = joinNodes ($2, NULL, NULL, NULL, T_NOT);}
	| expression AND expression {$$ = joinNodes ($1, $3, NULL, NULL, T_AND);}
	| expression OR expression {$$ = joinNodes ($1, $3, NULL, NULL, T_OR);}
	| OPEN expression CLOSE {$$ = $2;}
	;

unitspec: unitspec atomspec {$$ = joinNodes ($1, $2, NULL, NULL, T_OR);}
	| unitspec residuespec {$$ = joinNodes ($1, $2, NULL, NULL, T_OR);}
	| unitspec chainspec {$$ = joinNodes ($1, $2, NULL, NULL, T_OR);}
	| atomspec {$$ = $1;}
	| residuespec {$$ = $1;}
	| chainspec {$$ = $1;}
	;

atomspec: ATOM numlist {$2 = setTypeIntList ($2, T_ATOM); $$ = joinNodes (NULL, NULL, $2, NULL, T_GROUPING);}
	| ATOM namelist {$2 = setTypeNameList ($2, T_ATOM); /*printNameList ($2); printf("\n");*/ $$ = joinNodes (NULL, NULL, NULL, $2, T_GROUPING);}
	;

residuespec: RESIDUE numlist {$2 = setTypeIntList ($2, T_RESIDUE); $$ = joinNodes (NULL, NULL, $2, NULL, T_GROUPING);}
	| RESIDUE namelist {$2 = setTypeNameList ($2, T_RESIDUE); /*printNameList ($2); printf("\n");*/ $$ = joinNodes (NULL, NULL, NULL, $2, T_GROUPING);}
	;

chainspec: CHAIN chainlist {$2 = setTypeNameList ($2, T_CHAIN); /*printNameList ($2); printf("\n");*/ $$ = joinNodes (NULL, NULL, NULL, $2, T_GROUPING);}
	;

numlist: NUM {$$ = addIntIntList (NULL, $1);}
	| NUM RANGE NUM {$$ = addRangeIntList (NULL, $1, $3);}
	| numlist COMMA NUM {$$ = addIntIntList ($1, $3);}
	| numlist COMMA NUM RANGE NUM {$$ = addRangeIntList ($1, $3, $5);}
	;

namelist: NAME {/*printf("Name:%s\n", $1);*/ $$ = addNameNameList (NULL, $1);}
	| namelist COMMA NAME {$$ = addNameNameList ($1, $3);}
	;

chainlist: NAME {/*printf("Chainid:%s\n", $1);*/ $$ = addNameNameList (NULL, $1);}
  /*| NAME RANGE NAME {$$ = addRangeNameList (NULL, $1, $3);}*/
	| chainlist COMMA NAME {$$ = addNameNameList ($1, $3);}
	/*| chainlist COMMA NAME RANGE NAME {$$ = addRangeNameList ($1, $3, $5);}*/
	; 

%%

void yyerror (char const *s) {
	fprintf (stderr, "%s\n", s);
}
