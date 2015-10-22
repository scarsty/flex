FC=ifort
MPIF90=mpif90
FFLAGS=-I. -O2 -mcmodel=medium -shared-intel -w -mkl -fpp -u -extend-source
LINKFLAG=-L. -llapack -lblas
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)
TARGET=../ffff

.PHONY: all

all: $(TARGET)

mpi: FC=$(MPIF90) all

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(TARGET)

clean:
	rm *.o *.mod

%.o: %.f90
	$(FC) $(FFLAGS) -c $(@:.o=.f90)

# Dependencies of files
parameters2.o: \
    parameters2.f90 \
    constants.o
testmodule.o: \
    testmodule.f90 \
    constants.o \
    myfunctions.o \
    parameters2.o \
    fftw3.f03
myfunctions.o: \
    myfunctions.f90 \
    constants.o \
    parameters.o \
    parameters2.o \
    fftw3.f03
readin.o: \
    readin.f90 \
    constants.o \
    myfunctions.o \
    parameters.o \
    parameters2.o
main.o: \
    main.f90 \
    constants.o \
    myfunctions.o \
    parameters.o \
    parameters2.o
constants.o: \
    constants.f90
eliashberg.o: \
    eliashberg.f90 \
    constants.o \
    myfunctions.o \
    parameters.o \
    parameters2.o
parameters.o: \
    parameters.f90
