# makefile for the newton3 example 
# Here we use the macros CXX, CXXLAGS, and general rules for compiling 
# .cpp files.
# We also use the macro OBJECTS to not have to keep track of the 
# object files in more than one place.
# The PHONY target clean removes non-source fils.

AR = ar
.PHONY: clean d1arraylib d1arraylib_dbg d2arraylib d2arraylib_dbg

d1arraylib:
	g++ -O3 -c D1array.cpp
	$(AR) rcs libd1array.a D1array.o

d1arraylib_dbg: 
	g++ -g -DBZ_DEBUG -c D1array.cpp
	$(AR) rcs libd1array_dbg.a D1array.o

d2arraylib:
	g++ -O3 -c D2array.cpp
	$(AR) rcs libd2array.a D2array.o

d2arraylib_dbg:
	g++ -g -DBZ_DEBUG -c D2array.cpp
	$(AR) rcs libd2array_dbg.a D2array.o

clean:
	rm -f $(OBJECTS) $(EXEC)
