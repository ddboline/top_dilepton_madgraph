F77   = f77
FFLAGS= -O -ffixed-line-length-132
FFLAGS = -g -ffixed-line-length-132
LIBDIR = ../../lib/
BINDIR = ../../bin/
PROG   = madevent
LINKLIBS =  -L../../lib/ -ldhelas3 -ldsample -lgeneric -lpdf
LIBS = $(LIBDIR)libdhelas3.a $(LIBDIR)libdsample.a $(LIBDIR)libgeneric.a $(LIBDIR)libpdf.a


PROCESS= driver.o myamp.o auto_dsig.o genps.o matrix.o unwgt.o
SYMMETRY = myamp.o genps.o matrix.o 

%.o : %.f
	$(F77) $(FFLAGS) -c $<

$(PROG): $(PROCESS)  makefile $(LIBS)
	$(F77) $(FFLAGS) -o $(PROG) $(PROCESS) $(LINKLIBS)

symmetry.dat: mg.sym gensym
	./gensym

gensym: symmetry.o genps.inc $(SYMMETRY) configs.inc  
	$(F77) $(FFLAGS) -o gensym symmetry.o $(SYMMETRY) $(LINKLIBS)

driver.f: genps.inc 

driver.o: driver.f

genps.inc: nexternal.inc
	touch genps.inc

symmetry.o: genps.inc configs.inc
genps.o: genps.f genps.inc configs.inc 
cuts.o: cuts.f genps.inc pmass.inc
invarients.o: invarients.f genps.inc 
myamp.o: myamp.f props.inc genps.inc
unwgt.o: unwgt.f genps.inc symswap.inc
