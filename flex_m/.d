constants.o: constants.f90
eliashberg.o: eliashberg.f90 constants.o myfunctions.o parameters.o \
 parameters2.o
main.o: main.f90 constants.o myfunctions.o parameters.o \
 parameters2.o testmodule.o eliashberg.o readin.o
myfunctions.o: myfunctions.f90 constants.o parameters.o parameters2.o
parameters.o: parameters.f90
parameters2.o: parameters2.f90 constants.o
readin.o: readin.f90 constants.o parameters.o parameters2.o \
 myfunctions.o
testmodule.o: testmodule.f90 constants.o myfunctions.o parameters2.o
