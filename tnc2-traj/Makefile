.DEFAULT: .f .for .c .C .cpp .cc
.SUFFIXES: .f .for .c .C .cpp .cc

O = .

F77 = ifort
CC = icc
CCC = icpc 

LD = icpc

CFLAGS = -O3 -fPIC -align -Zp8 -axP -xP -unroll -Wno-deprecated -openmp -ip \
	-I/usr/include/python2.4 \
	-I/share/apps/python/2.4.3/lib64/python2.4/site-packages/numpy/core/include \
	-I/share/apps/netcdf/3.6.3/intel/include 

#	-I/home/sw77/QChem-Tinker/VMD-plugin 

LIBS = 	-L . \
	-L/share/apps/netcdf/3.6.3/intel/lib -lnetcdf_c++ -lnetcdf \
	-L/share/apps/intel/Compiler/11.1/046/lib/intel64 \
        -lifcore -lifport -lstdc++

EXENAME = QMMMMDtnc2.so
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
	cd $(O) ; $(CC)  $(CFLAGS) -c $<
$(O)/%.o: %.cc
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.cpp
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.C
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.for
	cd $(O) ; $(F77) $(FFLAGS) -c $<

clean:
	rm -f core *.o *.so
dat: 
	rm -f *.dat
backup:
	rm -f *~
clobber:
	rm -f $(EXENAME).exe
cleanall:
	rm -f *.o *.dat *~ *.exe *.exe.* $(EXENAME) *.pyc depend

depend :
	$(CCC) $(CFLAGS) -MM *.[cC] | perl dep.pl > $@


.PRECIOUS: %_wrap.C %.C %.c %.f %.h $(O)/%.o

include depend