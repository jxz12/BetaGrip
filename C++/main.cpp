#include <iostream>

#include "god.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "not enough arguments" << std::endl;
        return 0;
    }
    GodOfKeyboards god(argv[1]);
    auto result = god.BranchAndBound();
}