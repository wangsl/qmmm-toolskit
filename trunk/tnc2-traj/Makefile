
# $Id$ 

.DEFAULT: .f .for .c .C .cpp .cc
.SUFFIXES: .f .for .c .C .cpp .cc

O = .

SRC = /home/sw77/QChem4.0-Tinker4.2/qmmm-toolskit/tnc2-traj

O = /home/sw77/QChem4.0-Tinker4.2/qmmm-toolskit-build/tnc2-traj

F77 = ifort
CC = icc
CCC = icpc 

LD = icpc

CFLAGS = -O3 -fPIC -unroll -ip -axOP -xOP -openmp -vec-report -par-report -openmp-report -Wno-deprecated \
	-I/share/apps/python/2.7.2/intel/include/python2.7 \
	-I/share/apps/python/2.7.2/intel/lib/python2.7/site-packages/numpy/core/include \
	-I/share/apps/netcdf/3.6.3/intel/include 

FFLAGS = -O3 -fPIC -unroll -ip -axOP -xOP -openmp -vec-report -par-report -openmp-report

LIBS = 	-L/share/apps/python/2.7.2/intel/lib -lpython2.7 \
	-L/share/apps/netcdf/3.6.3/intel/lib -lnetcdf_c++ -lnetcdf \
	-L/share/apps/intel/Compiler/11.1/046/lib/intel64 \
        -lifcore -lifport -lstdc++

EXENAME = $(O)/QMMMMDtnc2.so
Link = $(LD) $(CFLAGS)

OBJS =  $(O)/extension.o  $(O)/readtnc2.o \
	$(O)/die.o  $(O)/out.o  $(O)/indent.o  $(O)/NetCDF.o \
	$(O)/NetCDFio.o \
	$(O)/ReadNetCDFtraj.o $(O)/ReadNetCDFtrajio.o \
	$(O)/utilC01.o

$(EXENAME) : $(OBJS) 
	$(Link) -shared -o $(EXENAME) $(OBJS) $(LIBS)  
	chmod -x $(EXENAME) 

$(O)/%.o: %.c
	cd $(O) ; $(CC)  $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.cc
	cd $(O) ; $(CCC) $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.cpp
	cd $(O) ; $(CCC) $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.C
	cd $(O) ; $(CCC) $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $(SRC)/$<
$(O)/%.o: %.for
	cd $(O) ; $(F77) $(FFLAGS) -c $(SRC)/$<

clean:
	rm -f $(O)/* *~ $(EXENAME) *.mod depend *.a *.pyc

depend :
	$(CCC) $(CFLAGS) -MM *.[cC] | perl dep.pl > $@


.PRECIOUS: %_wrap.C %.C %.c %.f %.h $(O)/%.o

#include depend	
