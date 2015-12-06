FC=ifort
MPIF90=mpif90
FFLAGS1=-I. -O3 -mcmodel=large -shared-intel -w -mkl -fpp -finline-functions
FFLAGS=$(FFLAGS1) -extend-source -no-inline-factor -qopenmp -parallel -xHost #-ipo
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
parameters.o: \
    parameters.f90
constants.o: \
    constants.f90
eliashberg.o: \
    eliashberg.f90 \
    constants.o \
    functions.o \
    parameters.o \
    parameters2.o
functions_test.o: \
    functions_test.f90 \
    functions_base.o \
    parameters.o \
    parameters2.o
functions_base.o: \
    functions_base.f90 \
    constants.o \
    parameters2.o \
    fftw3.f03
flex_2d.o: \
    flex_2d.f90 \
    constants.o \
    functions.o \
    parameters.o \
    parameters2.o
parameters2.o: \
    parameters2.f90 \
    constants.o
functions_init.o: \
    functions_init.f90 \
    functions_base.o \
    parameters.o \
    parameters2.o
readin.o: \
    readin.f90 \
    constants.o \
    functions.o \
    parameters.o \
    parameters2.o
functions_g_chi.o: \
    functions_g_chi.f90 \
    functions_base.o \
    parameters.o \
    parameters2.o
functions_mixer.o: \
    functions_mixer.f90 \
    functions_base.o \
    parameters.o \
    parameters2.o
functions.o: \
    functions.f90 \
    functions_base.o \
    functions_g_chi.o \
    functions_init.o \
    functions_mixer.o \
    functions_test.o

