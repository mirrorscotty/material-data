VPATH=../choi-okos/ matrix
CC=gcc
CFLAGS=-I../choi-okos -Imatrix -lm -ggdb
.PHONY: matrix

all: sens-analysis diff

composition.o: composition.c constants.h pasta.h choi-okos.h

thermal.o: thermal.c constants.h pasta.h choi-okos.h

fluid.o: fluid.c constants.h pasta.h isotherms.h

phase-change.o: phase-change.c constants.h pasta.h isotherms.h

gas.o: gas.c constants.h pasta.h

sensitivity.o: sensitivity.c pasta.h isotherms.h constants.h isotherms.h sensitivity.h

choi-okos.o: choi-okos.c choi-okos.h

isotherms.o: isotherms.c isotherms.h

sens-analysis: choi-okos.o isotherms.o composition.o thermal.o fluid.o phase-change.o gas.o sensitivity.o
	$(CC) -o sens-analysis choi-okos.o isotherms.o composition.o thermal.o fluid.o phase-change.o gas.o sensitivity.o $(CFLAGS)

matrix:
	$(MAKE) -C matrix

diffusivity.o: diffusivity.c diffusivity.h isotherms.h constants.h

binding.o: isotherms.h binding.c

diff-test.o: isotherms.h diffusivity.h matrix.h choi-okos.h diff-test.c

diff: matrix choi-okos.o diffusivity.o isotherms.o binding.o diff-test.o fluid.o composition.o gas.o
	$(CC) $(CFLAGS) -o diff -lm matrix/*.o choi-okos.o diffusivity.o isotherms.o binding.o fluid.o composition.o diff-test.o gas.o

clean:
	rm -rf *.o *.csv sensitivity sens-analysis diff


