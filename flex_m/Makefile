FC=ifort -I. -O2 -mcmodel=medium -shared-intel -w -mkl -fpp
MPIF90=mpif90 -I. -O2 -mcmodel=medium -shared-intel -w -mkl -fpp
FFLAGS=
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)
TARGET=ffff

.PHONY: all

all: $(TARGET)

mpi: FC=$(MPIF90)

$(TARGET): $(OBJECTS)
	$(FC) $(OBJECTS) -o $(TARGET)

clean:
	rm *.o *.mod

-include .d
