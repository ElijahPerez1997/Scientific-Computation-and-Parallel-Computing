#Makefile

results.csv: int.exe
	./int.exe > results.csv

int.exe: integral.o lglnodes.o
	gfortran integral.o lglnodes.o -o int.exe

integral.o: integral.f90
	gfortran -c integral.f90
lglnodes.o: lglnodes.f90
	gfortran -c lglnodes.f90