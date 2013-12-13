FC = mpif90.openmpi
FCFLAGS = -fbounds-check -ffree-form -pedantic-errors -Wall -Wconversion -Wunderflow -Wextra -Wunreachable-code -ffpe-trap=invalid,zero,overflow,underflow -g3 -fbounds-check -O0 -fstack-protector-all
#~ FCFLAGS = -O3 -std=f2003
FCFLAGS+= -std=f2003
PROGRAMS = schwarz_additif 
SRC= $(wildcard *.f90)
OBJ= $(SRC:.f90=.o)

all: $(PROGRAMS)

schwarz_additif: $(OBJ) utils.mod
	$(FC) $(FCFLAGS) -o $@ $(OBJ) $(LDFLAGS) >> make.log 2>&1

utils.mod : utilitaire.o utilitaire.f90
	$(FC) -std=f2003 -c $^

schw_add_MPI.o : schw_add_MPI.f90 utils.mod
	$(FC) $(FCFLAGS) -c $<

utilitaire.o : utilitaire.f90
	$(FC) -std=f2003 -c $<

second_membre.o : second_membre.f90
	$(FC) -std=f2003 -c $<

matvec.o : matvec.f90
	$(FC) -std=f2003 -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< >> make.log 2>&1

.PHONY: clean mrproper

clean:
	rm -f *.o *.mod *.MOD make.log

mrproper: clean
	rm -f *~ $(PROGRAMS)

