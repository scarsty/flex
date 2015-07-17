constants.o constants.mod: constants.f90
    $(FC) -c $(@:.o=.f90)
eliashberg.o: eliashberg.f90 constants.mod myfunctions.mod parameters.mod \
 parameters2.mod
    $(FC) -c $(@:.o=.f90)
main.o: main.f90 constants.mod myfunctions.mod parameters.mod \
 parameters2.mod
    $(FC) -c $(@:.o=.f90)
myfunctions.o myfunctions.mod: myfunctions.f90 constants.mod
    $(FC) -c $(@:.o=.f90)
parameters.o parameters.mod: parameters.f90
    $(FC) -c $(@:.o=.f90)
parameters2.o parameters2.mod: parameters2.f90 constants.mod
    $(FC) -c $(@:.o=.f90)
readin.o: readin.f90 constants.mod parameters.mod parameters2.mod \
 myfunctions.mod
    $(FC) -c $(@:.o=.f90)
