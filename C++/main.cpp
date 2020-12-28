#include <iostream>

#include "grip.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "not enough arguments" << std::endl;
        return 0;
    }
    BetaGrip grip(argv[1]);
    // auto result = grip.BranchAndBound();
    auto result = grip.SimulatedAnnealing(10000, .05, 1e3, .99);
}