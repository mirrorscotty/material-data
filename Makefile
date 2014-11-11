VPATH=composition diffusivity isotherms mechanical pasta test
CC=gcc
CFLAGS=-I. -Icomposition -Idiffusivity -Iisotherms -Imatrix -Imechanical -Ipasta -ggdb -O0
LDFLAGS=-lm

all: sens-analysis diff material-data.a pc_test

composition.o: composition.c constants.h pasta.h choi-okos.h
thermal.o: thermal.c constants.h pasta.h choi-okos.h
fluid.o: fluid.c constants.h pasta.h isotherms.h
phase-change.o: phase-change.c constants.h pasta.h isotherms.h
gas.o: gas.c constants.h pasta.h
sensitivity.o: sensitivity.c pasta.h isotherms.h constants.h isotherms.h

mechanical.o: mechanical.h
burgers.o: mechanical.h

diffusivity.o: diffusivity.h isotherms.h constants.h
gas-diff.o: diffusivity.h isotherms.h constants.h
capillary.o: diffusivity.h isotherms.h constants.h
binding.o: isotherms.h binding.c

oswin.o: isotherms.h
gab.o: isotherms.h
henderson.o: isotherms.h

choi-okos.o: choi-okos.h

diff-test.o: isotherms.h diffusivity.h matrix.a choi-okos.h diff-test.c
pc_test.o: isotherms.h diffusivity.h matrix.a choi-okos.h pc_test.c

material-data.a: composition.o thermal.o fluid.o phase-change.o gas.o choi-okos.o oswin.o gab.o henderson.o diffusivity.o capillary.o gas-diff.o binding.o mechanical.o burgers.o
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


