# Compiler
CXX = g++
MPICXX = mpic++

# Flags
CXXFLAGS_RELEASE = -O2 -Wall -I/usr/include/eigen3
CXXFLAGS_DEBUG = -g -O0 -DDEBUG -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -I/usr/include/eigen3
OPENMPFLAGS = -fopenmp
MPIFLAGS =

# Build type (default is Release)
ifeq ($(DEBUG),yes)
    CXXFLAGS = $(CXXFLAGS_DEBUG)
else
    CXXFLAGS = $(CXXFLAGS_RELEASE)
endif

# Targets
TARGETS = mi_serial mi_openmp mi_mpi

# Rules
all: $(TARGETS)

mi_serial: mi_serial.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

mi_openmp: mi_openmp.cpp
	$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -o $@ $<

mi_mpi: mi_mpi.cpp
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)
