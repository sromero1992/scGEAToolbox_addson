from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        name="mutual_information_serial",
        sources=["mutual_information_serial.pyx"],
        include_dirs=[numpy.get_include()]  # Add NumPy include path here
    )
]

setup(
    name="mutual_information_serial",
    ext_modules=cythonize(ext_modules),
)
