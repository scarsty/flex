constants.o : \
  constants.f90
	$(FC) -c $(@:.o=.f90)
main.o : \
  main.f90 constants.mod parameters2.F90 parameters.F90
	$(FC) -c $(@:.o=.f90)
mpi_functions.o : \
  mpi_functions.f90
	$(FC) -c $(@:.o=.f90)
readin.o : \
  readin.f90 mpi_functions.mod \
  constants.mod parameters.F90
	$(FC) -c $(@:.o=.f90)

