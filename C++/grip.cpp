#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cctype>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>

#include "grip.hpp"
#include "randomkit.h"

#define VERBOSE 1

inline unsigned BetaGrip::Two2OneD(unsigned row, unsigned col) {
    return row*CIRCUMF + col;
}

BetaGrip::BetaGrip(const std::string& textPath) {
    // read in text and count cooccurances between letters
    auto infile = std::ifstream(textPath);
    char c;
    while (infile.get(c)) {
        c = tolower(c);
        if (char2freq.find(c) == char2freq.end()) {
            char2freq[c] = 1;
        } else {
            char2freq[c] += 1;
        }
    }
    std::vector<char> idx2charAll;
    for (auto const& pair: char2freq) {
        idx2charAll.push_back(pair.first);
    }
    std::sort(idx2charAll.begin(), idx2charAll.end(), [this](
        char const& a, char const& b) {
            return char2freq[a] > char2freq[b];
        }
    );

    // will break if CIRCUMF > |set(letters)|
    for (unsigned i=0; i<CIRCUMF; i++) {
        char c = idx2charAll[i];
        idx2char[i] = c;
        char2idx[c] = i;
#if VERBOSE
        if (c == '\n') {
            std::cout << "\\n" << " " << char2freq[c] << std::endl;
        } else if (c == ' ') {
            std::cout << "\\s" << " " << char2freq[c] << std::endl;
        } else {
            std::cout << c << " " << char2freq[c] << std::endl;
        }
#endif
    }
    for (unsigned i=CIRCUMF; i<idx2charAll.size(); i++) {
        char c = idx2charAll[i];
        char2freq.erase(c);
    }


    // will break if '\n' not in letters
    std::fill(freqMatrix.begin(), freqMatrix.end(), 0);
    infile.clear();
    infile.seekg(0);
    char prev;
    infile.get(prev);
    prev = tolower(prev);
    auto end = char2idx.end();
    while (infile.get(c)) {
        c = tolower(c);
        if (char2idx.find(prev) != end && char2idx.find(c) != end) {
            auto row = char2idx[prev];
            auto col = char2idx[c];
            auto oneD = Two2OneD(row, col);
            freqMatrix[oneD] += 1;
        }
        prev = c;
    }
#if VERBOSE
    for (unsigned i=0; i<CIRCUMF; i++) {
        char c = idx2char[i];
        if (c == '\n') {
            std::cout << "\\n" << ": " << idx2char[i];
        } else if (c == ' ') {
            std::cout << "\\s" << ": " << char2freq[i];
        } else {
            std::cout << c << ":  " << char2freq[i];
        }
        for (unsigned j=0; j<CIRCUMF; j++) {
            auto k = Two2OneD(i, j);
            std::cout << freqMatrix[k] << '\t';
        }
        std::cout << std::endl;
    }
#endif
}

// an exhaustive search of all permutations
std::array<char, CIRCUMF> BetaGrip::BranchAndBound() {
    auto used = std::array<bool, CIRCUMF>();
    auto order = std::array<unsigned, CIRCUMF>();
    auto result = std::array<char, CIRCUMF>();
    unsigned best = -1;  // underflows to max int value
    
    std::function<void(unsigned, unsigned)> DFS;
    DFS = [&](unsigned pos, unsigned cost) {
        if (pos >= CIRCUMF) {
            for (unsigned i=0; i<CIRCUMF; i++) {
                result[i] = idx2char[order[i]];
                std::cout << (result[i]=='\n'? '<' : result[i]);
            }
            std::cout << ": " << cost << std::endl;
            best = cost;
            return;
        }
        for (unsigned i=0; i<CIRCUMF; i++) {
            if (!used[i]) {
                unsigned costNew = cost;
                for (unsigned pos2=0; pos2<pos; pos2++) {
                    auto dist = pos - pos2;
                    dist = std::min(dist, CIRCUMF-dist);  // clock/anticlockwise

                    auto j = order[pos2];
                    auto k = Two2OneD(i,j);
                    auto l = Two2OneD(j,i);
                    costNew += dist * (freqMatrix[k] + freqMatrix[l]);
                }
                if (costNew >= best) {
                    continue;  // quit early if this branch is a dead end
                }
                order[pos] = i;
                used[i] = true;
                DFS(pos+1, costNew);
                used[i] = false;
            }
        }
    };
    // fix first character
    order[0] = 0;
    result[0] = idx2char[0];
    used[0] = true;
    for (unsigned i=1; i<CIRCUMF; i++) {
        used[i] = false;
    }
    DFS(1,0);
    return result;
}

unsigned BetaGrip::EvalCost(std::array<unsigned, CIRCUMF> const& order) {
    unsigned cost = 0;
    for (unsigned pos=0; pos<CIRCUMF; pos++) {
        for (unsigned pos2=0; pos2<pos; pos2++) {
            auto dist = pos - pos2;
            dist = std::min(dist, CIRCUMF-dist);  // clock/anticlockwise

            auto i = order[pos];
            auto j = order[pos2];
            auto k = Two2OneD(i,j);
            auto l = Two2OneD(j,i);
            cost += dist * (freqMatrix[k] + freqMatrix[l]);
        }
    }
    return cost;
}

// simulated annealing
std::array<char, CIRCUMF> BetaGrip::SimulatedAnnealing(
    int nIter, double mutationProb, double tempInit, double tempCool, unsigned long rseed
) {
    // initialise
    auto order = std::array<unsigned, CIRCUMF>();
    auto result = std::array<char, CIRCUMF>();
    for (unsigned i=0; i<CIRCUMF; i++) {
        order[i] = i;
    }
    auto cost = EvalCost(order);

    rk_state rstate;
    rk_seed(rseed, &rstate);

    float temp = tempInit;
    for (unsigned iter=0; iter<nIter; iter++) {
        auto orderNew = order;
        // swap with random probability
        // a la fisher-yates
        for (unsigned i=CIRCUMF-1; i>0; i--) {
            auto rand = rk_double(&rstate);
            if (rand <= mutationProb) {
                auto j = rk_interval(i, &rstate);
                unsigned temp = orderNew[i];
                orderNew[i] = orderNew[j];
                orderNew[j] = temp;
            }
        }
        // accept new ordering if better
        // or if random acceptance function satisfied
        unsigned costNew = EvalCost(orderNew);
        bool update = false;
        if (costNew < cost) {
            update = true;
        } else {
            double costDelta = costNew - cost; // no worry about unsigned as delta>=0
            auto rand = rk_double(&rstate);
            if (exp(-costDelta/temp) > rand) {
                update = true;
            }
        }
        // update ordering depending on above
        if (update) {
            order = orderNew;
            cost = costNew;
            for (unsigned i=0; i<CIRCUMF; i++) {
                result[i] = idx2char[order[i]];
                std::cout << (result[i]=='\n'? '<' : result[i]);
            }
            std::cout << ": " << cost << " t=" << temp << std::endl;
        }
        // decay temperature
        temp *= tempCool;
    }
    return result;
}