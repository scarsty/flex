constants.o: constants.f90
	$(FC) -c $(@:.o=.f90)
eliashberg.o: eliashberg.f90 constants.o myfunctions.o parameters.o \
 parameters2.o
	$(FC) -c $(@:.o=.f90)
main.o: main.f90 constants.o myfunctions.o parameters.o \
 parameters2.o testmodule.o
	$(FC) -c $(@:.o=.f90)
myfunctions.o: myfunctions.f90 constants.o parameters.o parameters2.o
	$(FC) -c $(@:.o=.f90)
parameters.o: parameters.f90
	$(FC) -c $(@:.o=.f90)
parameters2.o: parameters2.f90 constants.o
	$(FC) -c $(@:.o=.f90)
readin.o: readin.f90 constants.o parameters.o parameters2.o \
 myfunctions.o
	$(FC) -c $(@:.o=.f90)
testmodule.o: testmodule.f90 constants.o myfunctions.o parameters2.o
	$(FC) -c $(@:.o=.f90) 
