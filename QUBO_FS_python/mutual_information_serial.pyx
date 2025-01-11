import numpy as np
cimport numpy as np
from scipy.stats import entropy
from scipy.sparse import csr_matrix

# A helper function to compute mutual information between two vectors
def compute_pairwise_mi(object matrix, int i, int j, int nbins):
    cdef np.ndarray[double, ndim=1] vi, vj
    cdef np.ndarray[double, ndim=2] joint_counts
    cdef double joint_prob, h_xy, h_x, h_y
    
    # Ensure the matrix is explicitly in double format
    matrix = matrix.astype('double')

    vi = matrix.getrow(i).toarray().flatten()
    vj = matrix.getrow(j).toarray().flatten()

    # Print types for debugging
    print("Type of vi:", vi.dtype)
    print("Type of vj:", vj.dtype)

    # Compute the 2D histogram (joint probability distribution)
    joint_counts, _, _ = np.histogram2d(vi, vj, bins=nbins)

    if joint_counts.sum() == 0:
        return 0  # No mutual information if no overlap

    joint_prob = joint_counts / (joint_counts.sum() + 1e-8)
    marginal_i = joint_prob.sum(axis=1) + 1e-8
    marginal_j = joint_prob.sum(axis=0) + 1e-8

    # Entropy calculation
    h_xy = entropy(joint_prob.flatten().astype('double'), base=2)
    h_x = entropy(marginal_i.astype('double'), base=2)
    h_y = entropy(marginal_j.astype('double'), base=2)

    return float(h_x + h_y - h_xy)

def mutual_information_matrix_cython(object matrix, int nbins=20):
    """
    Computes the mutual information matrix serially, working directly with sparse matrices,
    and only computes the upper triangular part of the matrix.
    This function is written in Cython for performance.
    """
    cdef int n_features = matrix.shape[0]
    cdef np.ndarray[double, ndim=2] mi_matrix = np.zeros((n_features, n_features), dtype='double')  # Change to double

    # Use serial loop for the upper triangular matrix, skip i == j
    for i in range(n_features):
        for j in range(i + 1, n_features):  # Start from i + 1 to avoid i == j
            mi_matrix[i, j] = compute_pairwise_mi(matrix, i, j, nbins)
            mi_matrix[j, i] = mi_matrix[i, j]  # Exploit symmetry

    return mi_matrix
