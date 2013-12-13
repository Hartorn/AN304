FC = mpif90.openmpi
FCFLAGS = -fbounds-check -ffree-form -std=f95 -pedantic-errors -Wall -Wconversion -Wimplicit-interface -Wunderflow -Wextra -Wunreachable-code -ffpe-trap=invalid,zero,overflow,underflow -g3 -fbounds-check -O0 -fstack-protector-all 
#~ FCFLAGS = -O3
PROGRAMS = schwarz_additif 
SRC= $(wildcard *.f90)
OBJ= $(SRC:.f90=.o)

all: $(PROGRAMS)

schwarz_additif: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
#~ 
#~ %: %.o
	#~ $(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean mrproper

clean:
	rm -f *.o *.mod *.MOD

mrproper: clean
	rm -f *~ $(PROGRAMS)

