#include <iostream>
#include <ctime>

#include "peach.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "not enough arguments" << std::endl;
        return 0;
    }
    PeachWheel grip(argv[1]);
    // auto result = grip.BruteForce();
    // auto result = grip.SimulatedAnnealing(1e6, 1e6, .99999, 0);
    auto result = grip.GeneticEvolution(10, 16, 2, 4, 0);
}