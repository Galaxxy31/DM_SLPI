F90=gfortran
EXE=exe
FLAGS=-pedantic -fcheck=all -Wall

#make + supression des fichiers *.o et *.mod
all : $(EXE) clean

#make + supression des fichiers *.o et *.mod + nettoyage terminal + execution
exec : $(EXE) clean clear execute

$(EXE) : mod_donnees.o mod_algorithmes.o main.o
	$(F90) -o $(EXE) $^
mod_donnees.o : mod_donnees.f90
	$(F90) -c $(FLAGS) $^
mod_algorithmes.o : mod_algorithmes.f90
	$(F90) -c $(FLAGS) $^
main.o : main.f90
	$(F90) -c $(FLAGS) $^

clean :
	rm -f *.o *.mod *.f90~

clear :
	clear

execute :
	./exe
