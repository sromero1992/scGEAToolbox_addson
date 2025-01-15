# Import necessary Cython features for Eigen and standard libraries
from libcpp.eigen cimport SparseMatrix, MatrixXd, VectorXd
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.algorithm cimport sort
from libcpp.math cimport log2
from libcpp.stdexcept cimport invalid_argument
from libc.stdlib cimport malloc, free
import numpy as np

# Declare the functions from your C++ code
cdef extern from "mi_serial.h":
    # C++ function declarations
    cdef MatrixXd calculate_mutual_information(const SparseMatrix[double] &sparseMatrix, size_t nbins)
    cdef double compute_pairwise_mi(const vector[double] &vec1, const vector[double] &vec2, size_t nbins)

# Python wrapper function to call the C++ mutual information calculation
def calculate_mutual_information_wrapper(np.ndarray[np.float64_t] sparse_matrix, size_t nbins):
    cdef int num_rows = sparse_matrix.shape[0]
    cdef int num_cols = sparse_matrix.shape[1]

    # Convert the NumPy array to Eigen SparseMatrix
    cdef SparseMatrix[double] cpp_sparse_matrix(num_rows, num_cols)

    # Fill Eigen SparseMatrix with the data from NumPy array
    for i in range(num_rows):
        for j in range(num_cols):
            cpp_sparse_matrix.insert(i, j) = sparse_matrix[i, j]

    # Call the C++ calculate_mutual_information function
    cdef MatrixXd cpp_mi_matrix = calculate_mutual_information(cpp_sparse_matrix, nbins)

    # Convert the resulting Eigen MatrixXd back to a NumPy array
    cdef np.ndarray[np.float64_t] result = np.zeros((num_rows, num_rows), dtype=np.float64)

    for i in range(num_rows):
        for j in range(num_rows):
            result[i, j] = cpp_mi_matrix(i, j)

    return result

