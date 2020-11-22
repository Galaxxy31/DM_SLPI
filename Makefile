F90=gfortran
EXE=exe
FLAGS=-pedantic -fcheck=all -Wall

1 : clear $(EXE) clean execute

$(EXE) : mod_donnees.o mod_algorithmes.o main.o
	$(F90) -o $(EXE) $^
mod_donnees.o : mod_donnees.f90
	$(F90) -c $(FLAGS) $^
mod_algorithmes.o : mod_algorithmes.f90
	$(F90) -c $(FLAGS) $^
main.o : main.f90
	$(F90) -c $(FLAGS) $^

clean :
	rm -f *.o *.mod

clear :
	clear

execute :
	./exe
