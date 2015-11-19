CC=gcc
CFLAGS=-I. -Icomposition -Idiffusivity -Iisotherms -Iglass-transition -Imath -Imatrix -Imechanical -Ipasta -ggdb -O0
LDFLAGS=-lm

SRC=$(wildcard test/*.c) $(LIBSRC)
LIBSRC=$(wildcard composition/*.c) \
       $(wildcard diffusivity/*.c) \
       $(wildcard glass-transition/*.c) \
       $(wildcard isotherms/*.c) \
       $(wildcard math/*.c) \
       $(wildcard mechanical/*.c) \
       $(wildcard pasta/*.c) \
       $(wildcard unifac/*.c)

all: sens-analysis diff material-data.a pc_test tg-test poisson-test

material-data.a: $(LIBSRC:.c=.o)
	ar -cvr $@ $?

matrix/matrix.a:
	$(MAKE) -C matrix

sens-analysis: test/pasta-sens.o material-data.a 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

diff: test/diff-test.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

pc_test: test/pc_test.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tg-test: test/tg-test.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

creep-test: test/creep-test.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

creep-testv2: test/creep-testv2.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

creep-tconst: test/creep-tconst.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

relax-test: test/relax-test.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

relax-exper: test/relax-exper.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

poisson-test: test/poisson-test.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

aw-calc: test/aw-calc.o material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

salt-calc: test/salt-calc.o material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

salt-gsp: test/salt-gsp.o material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

xdb-calc: test/xdb-calc.o material-data.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

aw-comp: test/aw-comp.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

aw-test: test/aw-test.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

aw-test-t: test/aw-test-t.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

plot-rozzi-t: plot-rozzi-t.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

plot-rozzi-xdb: test/plot-rozzi-xdb.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

isotherm-sensitivity: test/isotherm-sensitivity.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

diffusivity-sensitivity: test/diffusivity-sensitivity.o material-data.a matrix/matrix.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

ilt-test: test/ilt-test.o matrix/matrix.a material-data.a 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

co-test: test/co-test.o matrix/matrix.a material-data.a 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

unifac-data-test: test/unifac-data-test.o material-data.a matrix/matrix.a 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf sensitivity sens-analysis diff doc pc_test tg-test poisson-test relax-test creep-test
	rm -rf *.csv
	rm -rf *.a
	rm -rf $(SRC:.c=.d)
	rm -rf $(SRC:.c=.o)
	$(MAKE) -C matrix clean

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(CC) -MM $(CFLAGS) $*.c > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
          sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

-include $(SRC:.c=.d)

