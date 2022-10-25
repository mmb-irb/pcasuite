include ../config.conda

all: test_

# Linking targets
test_: test.o parser.tab.o lex.yy.o ast.o aux_parser.o ../utils.o

# Compilation targets
test.o: test.c y.tab.h aux_parser.h parser.h ast.h
aux_parser.o: aux_parser.c aux_parser.h ../pcz_io.h ../utils.h
ast.o: ast.c ast.h
parser.tab.o: parser.tab.c
y.tab.h: parser.y
	$(YACC) $(YFLAGS) parser.y
lex.yy.o: lex.yy.c y.tab.h
parser.tab.c: parser.y
	$(YACC) $(YFLAGS) parser.y
lex.yy.c : rules.f
	$(LEX) $(LFLAGS) rules.f

.PHONY: clean
clean:
	$(RM) -f *.o y.tab.h parser.tab.c parser.tab.bison lex.yy.c parser.output test
