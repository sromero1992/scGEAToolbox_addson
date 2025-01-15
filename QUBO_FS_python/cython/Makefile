CXX = g++
CFLAGS = -O2 -Wall -I/usr/include/eigen3 -I/usr/include/python3.11 -fPIC
LDFLAGS = -shared -lm

# Ensure that Python header files and libraries are correctly linked
PYTHON_INCLUDE = $(shell python3 -c 'from sysconfig import get_paths as gp; print(gp()["include"])')
PYTHON_LIB = $(shell python3 -c 'from sysconfig import get_paths as gp; print(gp()["stdlib"])')

all: mi_serial.so

mi_serial.so: mi_serial.o mi_serial_cython.o mmio.o
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lpython3.11

mi_serial.o: mi_serial.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

mi_serial_cython.o: mi_serial_cython.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

mmio.o: mmio.c
	cc -fPIC -c -o $@ $<

clean:
	rm -f *.o *.so

