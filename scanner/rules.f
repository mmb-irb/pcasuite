%option noyywrap
%option warn batch
%option nounput

%{
  #include <stdio.h>
	#include "ast.h"
	#include "y.tab.h"
%}

DIGIT	  [[:digit:]]
NAME	  [[:alnum:]]+

%%

[ \t\n]+	      /* Eat whitespace */
{DIGIT}+        yylval.num = atoi (yytext); return NUM;
\*?{NAME}\*?    yylval.name = yytext; return NAME;
-               return RANGE;
:               return RESIDUE;
@               return ATOM;
%               return CHAIN;
&               return AND;
\|              return OR;
~               return NOT;
\(              return OPEN;
\)              return CLOSE;
,               return COMMA;

%%

void processString (const char *string) {
	yy_scan_string (string);
}

void freeFlexResources (void) {
  yy_delete_buffer (YY_CURRENT_BUFFER);
}

/*
int main (int argc, char **argv) {
	yy_scan_string (":1-5&@CA|(~@H)");
	yylex();
	yy_delete_buffer (YY_CURRENT_BUFFER);
}
*/
