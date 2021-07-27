#Makefile for MPI version

FC = mpif90
OBJECTS = f2_fun.o f_mod.o laplac.o u_exact.o hw7_mpi.o

: hw7_mpi.exe
	./hw7_mpi.exe

hw7_mpi.exe: $(OBJECTS)
	$(FC) $(OBJECTS) -o hw7_mpi.exe

%.o : %.f90
	$(FC) -c $<

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -f hw4.exe
