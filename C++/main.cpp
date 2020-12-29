#include <iostream>
#include <ctime>

#include "peach.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "not enough arguments" << std::endl;
        return 0;
    }
    PeachWheel wheel(argv[1]);
    // auto result = wheel.BruteForce();
    auto result = wheel.RandomSample(1e7, 0);
    // auto result = wheel.SimulatedAnnealing(1e6, 1e6, .99999, 0);
    // auto result = wheel.GeneticEvolution(5000, 500, 5, 200, .5, 0);
}