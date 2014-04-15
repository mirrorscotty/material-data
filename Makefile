VPATH=choi-okos matrix pasta
CC=gcc
CFLAGS=-Ichoi-okos -Imatrix -ggdb
LDFLAGS=-lm

all: sens-analysis diff material-data.a

composition.o: composition.c constants.h pasta.h choi-okos.h

thermal.o: thermal.c constants.h pasta.h choi-okos.h

fluid.o: fluid.c constants.h pasta.h isotherms.h

phase-change.o: phase-change.c constants.h pasta.h isotherms.h

gas.o: gas.c constants.h pasta.h

sensitivity.o: sensitivity.c pasta.h isotherms.h constants.h isotherms.h

choi-okos.o: choi-okos.c choi-okos.h

isotherms.o: isotherms.c isotherms.h

diffusivity.o: diffusivity.c diffusivity.h isotherms.h constants.h

binding.o: isotherms.h binding.c

diff-test.o: isotherms.h diffusivity.h matrix.h choi-okos.h diff-test.c

material-data.a: composition.o thermal.o fluid.o phase-change.o gas.o choi-okos.o isotherms.o diffusivity.o binding.o
	ar -cvq $@ $?

matrix.a:
	$(MAKE) -C matrix
	cp matrix/matrix.a .

sens-analysis: sensitivity.o material-data.a 
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

diff: diff-test.o material-data.a matrix.a
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

clean:
	rm -rf *.o *.a *.csv sensitivity sens-analysis diff
	$(MAKE) -C matrix clean


