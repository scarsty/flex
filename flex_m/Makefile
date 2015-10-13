FC=ifort
MPIF90=mpif90
FFLAGS=-I. -O2 -mcmodel=medium -shared-intel -w -mkl -fpp -ffree-line-length-none
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

-include .d
