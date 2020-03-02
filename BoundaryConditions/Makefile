#
# Makefile for timing examples
#

CXX = g++
LD  = g++
SILO_LIB = /projects/mawa7160/APPM_HPC/Software/silo/lib
SILO_INC = /projects/mawa7160/APPM_HPC/Software/silo/include
DARRAY_LIB = ../darray_d
DARRAY_INC = ../darray_d
LDFLAGS = -L$(SILO_LIB) -lsilo -L$(DARRAY_LIB) -ldarray
CXXINCLUDE = -I$(DARRAY_INC) -I$(SILO_INC)
# CXXFLAGS = -g -O0 -pg -DBZ_DEBUG 
CXXFLAGS = -O3

OBJECTS = main.o BoundaryCondition.o
EXEC = ex.x	

.PHONY: clean

exe: $(EXEC)


$(EXEC): $(OBJECTS)
	$(LD) $(OBJECTS) -o $(EXEC) $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $<

clean:
	rm -f $(OBJECTS) $(EXEC) 
