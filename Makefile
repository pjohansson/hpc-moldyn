CC=g++

SRCDIR=src
TRGDIR=target
TESTDIR=tests

INCLUDEDIR=.

OPTS=-Wall -std=c++14 -I ${INCLUDEDIR}

main: conf.o
	${CC} ${OPTS} $^ ${SRCDIR}/main.cpp -o ${TRGDIR}/md

conf.o: ${SRCDIR}/conf.cpp
	${CC} -c ${OPTS} $< -o $@

tests: ${TESTDIR}/test_system_conf.cpp conf.o
	${CC} ${OPTS} $^ -o ${TRGDIR}/${TESTDIR}/test_system_init
