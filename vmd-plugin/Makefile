
# $Id: Makefile,v 1.1 2005/07/18 20:32:28 wangsl Exp wangsl $

.DEFAULT: .F .for .c .C .cpp .cc .f90
.SUFFIXES: .F .for .c .C .cpp .cc .f90

O = .

F77 = gfortran
CC = gcc
CCC = g++

FPP = fpp

CFLAGS = -O2 -fPIC -Wall -Wno-deprecated \
	-I/home/sw77/local/GNU/NetCDF-3.6.3/include

FFLAGS = -O2 -fPIC -Wall

Link = $(CCC) $(CFLAGS) -i-dynamic

LIBS = -L/home/sw77/local/GNU/NetCDF-3.6.3/lib \
	-lnetcdf_c++ -lnetcdf \
	-lstdc++ 

EXENAME = t1

SUBOBJS = $(O)/die.o  $(O)/out.o  $(O)/indent.o  $(O)/NetCDF.o \
        $(O)/NetCDFio.o \
	$(O)/ReadNetCDFtraj.o $(O)/ReadNetCDFtrajio.o \
	$(O)/utilC01.o

OBJS =  $(O)/main.o  $(SUBOBJS)

$(EXENAME) : $(OBJS) 
	ar -r libmyNetCDF.a $(SUBOBJS)
	 ranlib libmyNetCDF.a
	$(Link) -o $(EXENAME) $(OBJS) $(LIBS)

$(O)/%.o: %.c
	cd $(O) ; $(CC) $(CFLAGS) -c $<
$(O)/%.o: %.cc
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.cpp
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.C
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.F
	cd $(O); $(FPP) $*.F $*.f; $(F77) $(FFLAGS) -c $*.f; rm $*.f
$(O)/%.o: %.for
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.f90
	cd $(O) ; $(F90) $(FFLAGS) -c $<

cleanall:
	rm -f *.o *~ $(EXENAME) *.il *.mod depend *.a

depend :
	$(CCC) $(CFLAGS) -MM *.[cC] | perl dep.pl > $@


.PRECIOUS: %io.C %.C %.c %.f %.h $(O)/%.o

include depend
