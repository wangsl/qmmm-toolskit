
# The integer number in this Makefile is I4, 
# please uncomment this line in WHAM.py
# from numpy import int32 as intpy

.DEFAULT: .f .for .c .C .cpp .cc
.SUFFIXES: .f .for .c .C .cpp .cc

O = .

F77 = ifort
CC = icc
CCC = icpc 

CFLAGS = -O3 -fPIC -align -Zp8 -axP  -unroll -Wno-deprecated -openmp -ip -vec-report -par-report -openmp-report \
	-I/usr/include/python2.4 \
	-I/share/apps/python/2.4.3/lib64/python2.4/site-packages/numpy/core/include

FFLAGS = -O3 -fPIC -WB -align -Zp8 -axP -unroll -openmp -ip -vec-report -par-report -openmp-report  #-I8

LIBS = -L/share/apps/intel/Compiler/11.1/046/lib/intel64 \
	-limf -lifcore -lguide -lifport -lstdc++

Link = $(CC) $(CFLAGS)

EXENAME = whamit.so  

OBJS =  $(O)/wham.o  $(O)/whamit.o  $(O)/cputime.o  $(O)/header.o

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
	rm -f *.o *.dat *~ *.exe *.exe.* $(EXENAME) *.pyc

.PRECIOUS: %_wrap.C %.C %.c %.f %.h $(O)/%.o
