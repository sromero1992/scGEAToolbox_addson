#!/bin/bash

# Compiler
CXX=g++
MPICXX=mpic++

# Flags
CXXFLAGS="-O2 -Wall"
OPENMPFLAGS="-fopenmp"
MPIFLAGS=""

# Source Files
SERIAL_SRC="mi_serial.cpp"
OPENMP_SRC="mi_openmp.cpp"
MPI_SRC="mi_mpi.cpp"

# Targets
SERIAL_TARGET="mi_serial"
OPENMP_TARGET="mi_openmp"
MPI_TARGET="mi_mpi"

# Compile serial version
$CXX $CXXFLAGS -o $SERIAL_TARGET $SERIAL_SRC
echo "Compiled $SERIAL_TARGET"

# Compile OpenMP version
$CXX $CXXFLAGS $OPENMPFLAGS -o $OPENMP_TARGET $OPENMP_SRC
echo "Compiled $OPENMP_TARGET"

# Compile MPI version
$MPICXX $CXXFLAGS $MPIFLAGS -o $MPI_TARGET $MPI_SRC
echo "Compiled $MPI_TARGET"

