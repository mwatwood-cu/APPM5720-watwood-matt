# makefile for the Wave Equation 
# # Here we use the macros CXX, CXXLAGS, and general rules for compiling 
# # .cpp files.
# # We also use the macro OBJECTS to not have to keep track of the 
# # object files in more than one place.
# # The PHONY target clean removes non-source fils.

CXX = g++
LD  = g++

CXXINCLUDE = -I../darray
CXXFLAGS = -O0 -DBZ_DEBUG
EXEC = ./a.o

OBJECTS = main.o WaveSimulation.o
.PHONY: clean

$(EXEC): link $(OBJECTS)
	$(LD) $(OBJECTS) $(LDFLAGS) -o $(EXEC)

link:
ifneq (,$(findstring DBZ_DEBUG,$(CXXFLAGS)))
	$(eval LDFLAGS := -L../darray -ld1array_dbg -ld2array_dbg)
	echo $(LDFLAGS)
else
	$(eval LDFLAGS := -L../darray -ld1array -ld2array)
	echo $(LDFLAGS)
endif

run_it: $(EXEC)
	$(EXEC)   

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $<

clean:
	rm -f $(OBJECTS) $(EXEC)
