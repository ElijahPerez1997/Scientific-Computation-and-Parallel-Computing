#Makefile

FC = gfortran
FFLAGS = -O3
OBJECTS = differentiate.o printdble2d.o xycoord.o trap2d.o homework4.o

: hw4.exe
	./hw4.exe

hw4.exe: $(OBJECTS)
	$(FC) $(OBJECTS) -o hw4.exe

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -f hw4.exe
