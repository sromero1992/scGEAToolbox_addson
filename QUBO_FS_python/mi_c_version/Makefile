# Compiler
CXX = g++
MPICXX = mpic++

# Flags
CXXFLAGS_RELEASE = -O2 -Wall -I/usr/include/eigen3 -lm
CXXFLAGS_DEBUG = -g -O0 -DDEBUG -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -I/usr/include/eigen3 -lm

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

mi_serial: mi_serial.cpp mmio.c mmio.h
	$(CXX) $(CXXFLAGS) -o $@ mi_serial.cpp mmio.c

mi_openmp: mi_openmp.cpp mmio.c mmio.h
	$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -o $@ mi_openmp.cpp mmio.c

mi_mpi: mi_mpi.cpp mmio.c mmio.h
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ mi_mpi.cpp mmio.c

clean:
	rm -f $(TARGETS)

