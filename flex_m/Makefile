FC=ifort
MPIF90=mpif90
FFLAGS=-I. -O3 -mcmodel=large -shared-intel -w -mkl -fpp -extend-source -finline-functions -no-inline-factor -qopenmp
MPIFFLAGS=$(FFLAGS) -DUSE_MPI
#LDFLAGS=-L. -llapack -lblas
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)
TARGET=../ffff

.PHONY: all

all: $(TARGET)

mpi: FC=$(MPIF90) all

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

clean:
	rm *.o *.mod

%.o: %.f90
	$(FC) $(FFLAGS) -c $(@:.o=.f90)

# Dependencies of files
eliashberg.o: \
    eliashberg.f90 \
    constants.o \
    functions.o \
    parameters.o \
    parameters2.o
flex_2d.o: \
    flex_2d.f90 \
    constants.o \
    functions.o \
    parameters.o \
    parameters2.o
readin.o: \
    readin.f90 \
    constants.o \
    functions.o \
    parameters.o \
    parameters2.o
functions.o: \
    functions.f90 \
    constants.o \
    parameters.o \
    parameters2.o \
    fftw3.f03
parameters2.o: \
    parameters2.f90 \
    constants.o
testmodule.o: \
    testmodule.f90 \
    constants.o \
    functions.o \
    parameters2.o \
    fftw3.f03
parameters.o: \
    parameters.f90
constants.o: \
    constants.f90
