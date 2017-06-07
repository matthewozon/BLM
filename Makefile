F90 = gfortran
FFLAGS = -O0 -g -ffpe-trap=invalid,zero,overflow -ffree-line-length-none

OBJS = parameters_mod.o time_mod.o grid_mod.o meteorology_mod.o main.o
EXE = main.exe

all: $(EXE)

$(EXE): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

clean:
	@rm -v -f *.o *.mod main.exe
