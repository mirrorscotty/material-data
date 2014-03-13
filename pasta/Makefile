VPATH=../choi-okos/ ../diffusivity
CC=gcc
CFLAGS=-I../choi-okos -I../diffusivity

all: sens-analysis

composition.o: composition.c constants.h pasta.h choi-okos.h

thermal.o: thermal.c constants.h pasta.h choi-okos.h

fluid.o: fluid.c constants.h pasta.h isotherms.h diff_data.h

phase-change.o: phase-change.c constants.h pasta.h isotherms.h diff_data.h

gas.o: gas.c constants.h pasta.h

sensitivity.o: sensitivity.c pasta.h isotherms.h constants.h isotherms.h diff_data.h sensitivity.h

choi-okos.o: choi-okos.c choi-okos.h

isotherms.o: isotherms.c isotherms.h

diff_data.o: diff_data.c diff_data.h

sens-analysis: diff_data.o choi-okos.o isotherms.o composition.o thermal.o fluid.o phase-change.o gas.o sensitivity.o
	gcc -o sens-analysis -lm diff_data.o choi-okos.o isotherms.o composition.o thermal.o fluid.o phase-change.o gas.o sensitivity.o

clean:
	rm -rf *.o sensitivity sens-analysis


