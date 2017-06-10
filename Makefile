F90 = gfortran-5
FFLAGS = -O0 -g -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -std=legacy

OBJS = parameters_mod.o time_mod.o grid_mod.o chemf.o opkdmain.o opkda1.o opkda2.o meteorology_mod.o main.o
EXE = blm

all: $(EXE)

#general rules
$(EXE): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $^

#rule for the f90 files
%.o: %.f90
	$(F90) $(FFLAGS) -c $<


# special rules for the chemistry
chemf.o: chemf.f90 opkdmain.o opkda1.o opkda2.o
	$(F90) $(FFLAGS) -c chemf.f90

opkda1.o: opkda1.f
	$(F90) $(FFLAGS)  -w -c opkda1.f

opkda2.o: opkda2.f
	$(F90) $(FFLAGS)  -c opkda2.f

opkdmain.o: opkdmain.f
	$(F90) $(FFLAGS)  -c opkdmain.f

clean:
	@rm -v -f *.o *.mod blm
