#
#   Makefile
#
FC = gfortran
FFLAGS = -Og  --pedantic -ftrapv -fbacktrace -g -fcheck=all -fdefault-integer-8
LIBS = -lblas -llapack 
#LIBS = -L/opt/OpenBLAS/lib -lopenblas

MODS   = real_precision.o utils.o diaglib.o
OBJS   = main.o
#
all:    $(MODS) $(OBJS)
	$(FC) $(FFLAGS) -o main.exe $(OBJS) $(MODS) $(LIBS)
#
%.o: %.f
	$(FC) $(FFLAGS) -c $*.f
%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(MODS) $(OBJS) *.exe *.mod
