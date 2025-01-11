#include <iostream>
#include <vector>
#include <cmath>
#include <map>

double mutualInformation(const std::vector<int>& X, const std::vector<int>& Y) {
    if (X.size() != Y.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    // Use a map of pairs to store joint frequencies
    std::map<std::pair<int, int>, int> freqXY; 
    std::map<int, int> freqX, freqY;

    // Calculate frequencies
    for (size_t i = 0; i < X.size(); ++i) {
        freqX[X[i]]++;
        freqY[Y[i]]++;
        freqXY[{X[i], Y[i]}]++; // Use {X[i], Y[i]} to create a pair
    }

    double mi = 0.0;
    for (const auto& pair : freqXY) {
        int x = pair.first.first;
        int y = pair.first.second;
        double px = static_cast<double>(freqX[x]) / X.size();
        double py = static_cast<double>(freqY[y]) / Y.size();
        double pxy = static_cast<double>(pair.second) / X.size();

        mi += pxy * log2(pxy / (px * py));
    }

    return mi;
}

int main() {
    std::vector<int> X = {1, 1, 0, 0, 1};
    std::vector<int> Y = {0, 1, 0, 1, 0};

    double mi = mutualInformation(X, Y);

    std::cout << "Mutual Information: " << mi << std::endl;

    return 0;
}