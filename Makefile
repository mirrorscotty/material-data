VPATH=choi-okos matrix pasta test pasta/isotherms
CC=gcc
CFLAGS=-Ichoi-okos -Imatrix -Ipasta -Ipasta/isotherms -ggdb
LDFLAGS=-lm

all: sens-analysis diff material-data.a pc_test

composition.o: composition.c constants.h pasta.h choi-okos.h
thermal.o: thermal.c constants.h pasta.h choi-okos.h
fluid.o: fluid.c constants.h pasta.h isotherms.h
phase-change.o: phase-change.c constants.h pasta.h isotherms.h
gas.o: gas.c constants.h pasta.h
sensitivity.o: sensitivity.c pasta.h isotherms.h constants.h isotherms.h
choi-okos.o: choi-okos.c choi-okos.h
mechanical.o: mechanical.h
burgers.o: mechanical.h
diffusivity.o: diffusivity.c diffusivity.h isotherms.h constants.h
binding.o: isotherms.h binding.c

oswin.o: isotherms.h
gab.o: isotherms.h
henderson.o: isotherms.h

diff-test.o: isotherms.h diffusivity.h matrix.h choi-okos.h diff-test.c
pc_test.o: isotherms.h diffusivity.h matrix.h choi-okos.h pc_test.c

material-data.a: composition.o thermal.o fluid.o phase-change.o gas.o choi-okos.o oswin.o gab.o henderson.o diffusivity.o binding.o mechanical.o burgers.o
	ar -cvr $@ $?

matrix.a:
	$(MAKE) -C matrix
	cp matrix/matrix.a .

sens-analysis: sensitivity.o material-data.a 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

diff: diff-test.o material-data.a matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

pc_test: pc_test.o material-data.a matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o *.a *.csv sensitivity sens-analysis diff doc pc_test
	$(MAKE) -C matrix clean


