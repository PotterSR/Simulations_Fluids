FC     = gfortran
FFLAGS = -O2 -Wall -fopenmp

OBJS = helpers.o params.o functions.o montecarlo.o

sim: $(OBJS)
	$(FC) $(FFLAGS) -o sim $(OBJS)

helpers.o: helpers.f90
	$(FC) $(FFLAGS) -c helpers.f90

params.o: params.f90 helpers.o
	$(FC) $(FFLAGS) -c params.f90

functions.o: functions.f90 helpers.o params.o
	$(FC) $(FFLAGS) -c functions.f90

montecarlo.o: montecarlo.f90 helpers.o params.o functions.o
	$(FC) $(FFLAGS) -c montecarlo.f90

clean:
	rm -f *.o *.mod sim