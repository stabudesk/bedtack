CC=gcc
CFLAGS=-O3
DBGCFLAGS=-g -Wall
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!

EXES=bedtack

# production binary
bedtack: bedtack.c
	${CC} ${CFLAGS} -o $@ $^

# testing mode binary
bedtack_t: bedtack.c
	${CC} ${DBGCFLAGS} -o $@ $^

# testing mode binary
bedtack_d: bedtack.c
	${CC} ${TDBGCFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
