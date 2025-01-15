from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        name="mi_serial_cython",  # Name of the compiled module
        sources=["mi_serial_cython.pyx", "mi_serial.cpp", "mmio.c"],
        language="c++",
        extra_compile_args=["-O2", "-fopenmp"],  # Optimize and use OpenMP if applicable
        include_dirs=[numpy.get_include(), "/usr/include/eigen3"],  # Add Eigen path
    )
]

setup(
    name="mi_serial_cython",
    ext_modules=cythonize(extensions),
    include_dirs=[numpy.get_include()],
)

