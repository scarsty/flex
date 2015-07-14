constants.o : \
  constants.f90
	$(FC) -c $(@:.o=.f90)
main.o : \
  main.f90 constants.o parameters2.F90 parameters.F90 myfunctions.o
	$(FC) -c $(@:.o=.f90)
myfunctions.o : \
  mpi_functions.f90
	$(FC) -c $(@:.o=.f90)
readin.o : \
  readin.f90 myfunctions.o \
  constants.o parameters.F90
	$(FC) -c $(@:.o=.f90)

