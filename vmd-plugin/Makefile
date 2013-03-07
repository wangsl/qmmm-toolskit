
# $Id$

.DEFAULT: .F .for .c .C .cpp .cc .f90
.SUFFIXES: .F .for .c .C .cpp .cc .f90

SRC = /home/sw77/QChem4.0-Tinker4.2/qmmm-toolskit/vmd-plugin

O = /home/sw77/QChem4.0-Tinker4.2/qmmm-toolskit-build/vmd-plugin

F77 = ifort
CC = icc
CCC = icpc

CFLAGS = -O3 -fPIC -unroll -ip -axOP -xOP -openmp -vec-report -par-report -openmp-report -Wno-deprecated \
	-I/home/sw77/QChem4.0-Tinker4.2/ext-libs/include 

FFLAGS = -O3 -fPIC -unroll -ip -axOP -xOP -openmp -vec-report -par-report -openmp-report 

Link = $(CCC) $(CFLAGS) -shared-intel

LIBS = -L/home/sw77/QChem4.0-Tinker4.2/ext-libs/lib \
	-lnetcdf_c++ -lnetcdf \
	-lstdc++ 

all: $(O)/t1 $(O)/TinkerNetCDFplugin.so  $(O)/TinkerNetCDF2plugin.so

EXENAME = $(O)/t1 

SUBOBJS = $(O)/die.o  $(O)/out.o  $(O)/indent.o  $(O)/NetCDF.o \
        $(O)/NetCDFio.o \
	$(O)/ReadNetCDFtraj.o $(O)/ReadNetCDFtrajio.o \
	$(O)/utilC01.o

OBJS =  $(O)/main.o  $(SUBOBJS)

$(EXENAME) : $(O) $(OBJS) 
	$(Link) -o $(EXENAME) $(OBJS) $(LIBS)

$(O)/%.o: %.c
	cd $(O) ; $(CC) $(CFLAGS) -c  $(SRC)/$<
$(O)/%.o: %.cc
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.cpp
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.C
	cd $(O) ; $(CCC) $(CFLAGS) -c  $(SRC)/$<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $<	
$(O)/%.o: %.for
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.f90
	cd $(O) ; $(F90) $(FFLAGS) -c $<

$(O) :
	mkdir -p $(O)

$(O)/%.so : %.C $(SUBOBJS) 
	cd $(O) ; $(Link) -shared -o $@  $(SRC)/$< $(SUBOBJS) $(LIBS)

clean:
	rm -f $(O)/* *~ $(EXENAME) *.mod depend *.a

depend :
	$(CCC) $(CFLAGS) -MM *.[cC] | perl dep.pl > $@

.PRECIOUS: %io.C %.C %.c %.f %.h $(O)/%.o


