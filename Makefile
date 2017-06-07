F90 = gfortran-5
FFLAGS = -O0 -g -ffpe-trap=invalid,zero,overflow -ffree-line-length-none

OBJS = parameters_mod.o time_mod.o grid_mod.o meteorology_mod.o main.o
EXE = blm

all: $(EXE)

$(EXE): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

clean:
	@rm -v -f *.o *.mod blm
