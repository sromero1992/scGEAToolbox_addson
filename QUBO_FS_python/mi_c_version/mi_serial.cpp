#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra> // For Matrix Market file support
#include <limits>
#include <stdexcept>

// Entropy calculation
double entropy(const std::vector<double>& prob) {
    if (prob.empty()) {
        std::cerr << "Error: Probability vector is empty!" << std::endl;
        return 0.0;
    }

    std::vector<double> prob_copy = prob; // Create a copy of the vector
    double total_prob = std::accumulate(prob_copy.begin(), prob_copy.end(), 0.0);
    if (total_prob > 0) {
        for (double& p : prob_copy) {
            p /= total_prob;  // Normalize to ensure sum is 1
        }
    }

    double result = 0.0;
    for (double p : prob_copy) {
        if (p > 0) {
            result -= p * std::log2(p);
        }
    }
    return result;
}

// Helper function for quantile binning
int compute_bin_index_quantile(double value, const std::vector<double>& bin_edges) {
    for (size_t i = 0; i < bin_edges.size() - 1; ++i) {
        if (value >= bin_edges[i] && value < bin_edges[i + 1]) {
            return static_cast<int>(i); // Cast to int for return type
        }
    }
    // Ensure values equal to the last edge are assigned to the last bin
    if (value == bin_edges.back()) {
        return static_cast<int>(bin_edges.size() - 2); // Last valid bin index
    }
    return static_cast<int>(bin_edges.size() - 1); // Return the last bin (in case of float precision issues)
}

double compute_pairwise_mi(const Eigen::SparseMatrix<double>& sparse_matrix, int i, int j, size_t nbins) {
    std::cout << "Computing pairwise MI for indices i = " << i << " and j = " << j << std::endl;

    // Convert sparse matrix row i and j to full vectors (including zeros)
    Eigen::VectorXd values_i(sparse_matrix.cols());
    Eigen::VectorXd values_j(sparse_matrix.cols());

    values_i.setZero();  // Initialize the vectors with zeroes
    values_j.setZero();

    bool has_nonzero_i = false, has_nonzero_j = false;

    // Populate the vectors with non-zero values from sparse matrix
    for (int k = 0; k < sparse_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(sparse_matrix, k); it; ++it) {
            if (it.row() == i) {
                values_i[it.col()] = it.value();  // Store non-zero value in dense vector
                has_nonzero_i = true;
            }
            if (it.row() == j) {
                values_j[it.col()] = it.value();  // Store non-zero value in dense vector
                has_nonzero_j = true;
            }
        }
    }

    // Print the values for debugging
    std::cout << "Values for i (" << i << "): ";
    for (int idx = 0; idx < values_i.size(); ++idx) {
        std::cout << values_i[idx] << " ";
    }
    std::cout << std::endl;

    std::cout << "Values for j (" << j << "): ";
    for (int idx = 0; idx < values_j.size(); ++idx) {
        std::cout << values_j[idx] << " ";
    }
    std::cout << std::endl;

    // Optionally, check if any values were populated
    if (!has_nonzero_i) {
        std::cout << "No non-zero values found for row " << i << std::endl;
    }
    if (!has_nonzero_j) {
        std::cout << "No non-zero values found for row " << j << std::endl;
    }

    // If either row i or row j is full of zeros, skip MI computation for this pair
    if (!has_nonzero_i && !has_nonzero_j) {
        std::cout << "The vectors are zero, returning MI = 0 for i = " << i << " and j = " << j << std::endl;
        return 0.0;
    }

    // Create quantile bin edges for i and j
    std::vector<double> bin_edges_i(nbins + 1), bin_edges_j(nbins + 1);

    // Sort the values and compute bin edges based on sorted data
    std::vector<double> sorted_i(values_i.data(), values_i.data() + values_i.size());
    std::vector<double> sorted_j(values_j.data(), values_j.data());
    std::sort(sorted_i.begin(), sorted_i.end());
    std::sort(sorted_j.begin(), sorted_j.end());

    // Compute the bin edges based on quantiles
    for (size_t k = 0; k <= nbins; ++k) {
        bin_edges_i[k] = sorted_i[static_cast<size_t>(k * sorted_i.size() / nbins)];
        bin_edges_j[k] = sorted_j[static_cast<size_t>(k * sorted_j.size() / nbins)];
    }
    std::cout << "Quantile bin edges for i and j computed." << std::endl;

    // Create a 2D histogram (joint probability distribution)
    std::vector<std::vector<int>> joint_counts(nbins, std::vector<int>(nbins, 0));

    std::cout << "Starting joint histogram computation." << std::endl;

    // Compute the joint histogram
    for (int k = 0; k < sparse_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(sparse_matrix, k); it; ++it) {
            if (it.row() == i || it.row() == j) {
                int bin_i = compute_bin_index_quantile(it.value(), bin_edges_i);
                int bin_j = compute_bin_index_quantile(it.value(), bin_edges_j);
                joint_counts[bin_i][bin_j]++;
            }
        }
    }

    std::cout << "Joint histogram computation completed." << std::endl;

    // Now calculate the mutual information from the joint histogram
    double result = 0.0;
    int total_count = 0;
    for (size_t bin_i = 0; bin_i < nbins; ++bin_i) {
        for (size_t bin_j = 0; bin_j < nbins; ++bin_j) {
            total_count += joint_counts[bin_i][bin_j];
        }
    }

    std::cout << "Total count in joint histogram: " << total_count << std::endl;

    // Calculate the mutual information
    for (size_t bin_i = 0; bin_i < nbins; ++bin_i) {
        for (size_t bin_j = 0; bin_j < nbins; ++bin_j) {
            double p_ij = static_cast<double>(joint_counts[bin_i][bin_j]) / total_count;
            double p_i = static_cast<double>(std::accumulate(joint_counts[bin_i].begin(), joint_counts[bin_i].end(), 0)) / total_count;
            double p_j = static_cast<double>(std::accumulate(joint_counts.begin(), joint_counts.end(), 0, 
                        [bin_j](int acc, const std::vector<int>& row) { return acc + row[bin_j]; })) / total_count;

            if (p_ij > 0) {
                result += p_ij * log(p_ij / (p_i * p_j));
            }
        }
    }

    return result;  // Return the calculated MI
}




// Load sparse matrix from Matrix Market file
Eigen::SparseMatrix<double> loadMatrixMarketFile(const std::string& filename) {
    std::cout << "Loading matrix from file: " << filename << std::endl;

    Eigen::SparseMatrix<double> matrix;
    if (!Eigen::loadMarket(matrix, filename)) {
        throw std::runtime_error("Failed to load the Matrix Market file: " + filename);
    }

    std::cout << "Matrix loaded with size: " << matrix.rows() << "x" << matrix.cols() << std::endl;

    // Debugging: Print a few matrix entries
    int max_entries_to_print = 10; // Limit the number of entries printed
    int count = 0;

    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
            if (count >= max_entries_to_print) {
                break; // Stop printing after max_entries_to_print
            }
            std::cout << "Matrix entry: (" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
            count++;
        }
        if (count >= max_entries_to_print) {
            break; // Stop iterating after max_entries_to_print
        }
    }

    return matrix;
}

// Compute error between two matrices
double computeError(const Eigen::MatrixXd& matrix1, const Eigen::MatrixXd& matrix2) {
    if (matrix1.rows() != matrix2.rows() || matrix1.cols() != matrix2.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for error computation.");
    }

    double error = 0.0;
    for (int i = 0; i < matrix1.rows(); ++i) {
        for (int j = 0; j < matrix1.cols(); ++j) {
            error += std::pow(matrix1(i, j) - matrix2(i, j), 2); // Mean Squared Error
        }
    }
    return std::sqrt(error / (matrix1.rows() * matrix1.cols())); // Root Mean Squared Error
}

// Main function
int main() {
    Eigen::SparseMatrix<double> sparseMatrix;
    Eigen::SparseMatrix<double> pythonMI;

    try {
        // Load sparse input matrix from Matrix Market file
        sparseMatrix = loadMatrixMarketFile("sparse_matrix.mtx");

        // Load Python-computed MI matrix
        pythonMI = loadMatrixMarketFile("mi_matrix.mtx");

        std::cout << "Sparse Matrix size: " << sparseMatrix.rows() << "x" << sparseMatrix.cols() << std::endl;
        std::cout << "Python MI Matrix size: " << pythonMI.rows() << "x" << pythonMI.cols() << std::endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    // Compute MI matrix in C++ using sparse matrix directly
    size_t nbins = 20; // Number of bins for discretization
    Eigen::MatrixXd cppMI(sparseMatrix.rows(), sparseMatrix.rows()); // Use cols for MI matrix

    for (int i = 0; i < sparseMatrix.rows(); ++i) {  // Iterate over columns
        for (int j = i; j < sparseMatrix.rows(); ++j) {  // Iterate over columns
            double mi_value = compute_pairwise_mi(sparseMatrix, i, j, nbins);
            cppMI(i, j) = mi_value;
            cppMI(j, i) = mi_value; // Since MI is symmetric
        }
    }

    // Calculate error between C++ and Python MI matrices
    try {
        double error = computeError(cppMI, pythonMI);
        std::cout << "Error between C++ and Python MI matrices: " << error << std::endl;
    } catch (const std::invalid_argument& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }

    return 0;
}