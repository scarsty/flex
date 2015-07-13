F90=ifort -I. -O2 -mcmodel=medium -shared-intel -w -mkl -fpp
MPIF90=mpif90 -I. -O2 -mcmodel=medium -shared-intel -w -mkl -fpp
FCFLAGS=
SOURCES=constants.f90 mpi_functios.f90 readin.f90 main.f90
OBJECTS=$(SOURCES:.f90=.o)
TARGET=ffff

.PHONY: all

all: $(TARGET)

$(OBJECTS): $(SOURCES)
	$(F90) -c $(SOURCES)

$(TARGET): $(OBJECTS)
	$(F90) $(OBJECTS) -o $(TARGET)

clean:
	rm *.o .d *.mod
