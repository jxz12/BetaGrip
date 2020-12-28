#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <list>
#include <array>
#include <string>
#include <algorithm>
#include <functional>
#include <cctype>
#include <cmath>

#include "grip.hpp"
#include "randomkit.h"

#define VERBOSE 1

inline uint BetaGrip::Two2OneD(uint const& row, uint const& col) {
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
    std::sort(idx2charAll.begin(), idx2charAll.end(), [this](char a, char b) {
        return char2freq[a] > char2freq[b];
    });

    // will break if CIRCUMF > |set(letters)|
    for (uint i=0; i<CIRCUMF; i++) {
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
    for (uint i=CIRCUMF; i<idx2charAll.size(); i++) {
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
    for (uint i=0; i<CIRCUMF; i++) {
        char c = idx2char[i];
        if (c == '\n') {
            std::cout << "\\n" << ": " << idx2char[i];
        } else if (c == ' ') {
            std::cout << "\\s" << ": " << char2freq[i];
        } else {
            std::cout << c << ":  " << char2freq[i];
        }
        for (uint j=0; j<CIRCUMF; j++) {
            auto k = Two2OneD(i, j);
            std::cout << freqMatrix[k] << '\t';
        }
        std::cout << std::endl;
    }
#endif
}

inline std::string BetaGrip::Idxs2String(std::array<uint, CIRCUMF> const& order) {
    auto str = std::string(CIRCUMF, ' ');
    for (uint i=0; i<CIRCUMF; i++) {
        char c = idx2char[order[i]];
        str[i] = c;
    }
#if VERBOSE
    for (uint i=0; i<CIRCUMF; i++) {
        std::cout << (str[i]=='\n'? '_':str[i]);
    }
    std::cout << std::endl;
#endif
    return str;
}

// an exhaustive search of all permutations
std::string BetaGrip::BruteForce() {
    auto used = std::array<bool, CIRCUMF>();
    std::fill(used.begin(), used.end(), false);

    auto order = std::array<uint, CIRCUMF>();
    auto result = std::string(CIRCUMF, ' ');
    uint best = -1;  // underflows to max int value
    
    std::function<void(uint, uint)> DFS;
    DFS = [&](uint pos, uint cost) {
        if (pos >= CIRCUMF) {
            result = Idxs2String(order);
            best = cost;
            return;
        }
        for (uint i=0; i<CIRCUMF; i++) {
            if (!used[i]) {
                uint costNew = cost;
                for (uint pos2=0; pos2<pos; pos2++) {
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
    used[0] = true;
    DFS(1,0);
    return result;
}

uint BetaGrip::EvalCost(std::array<uint, CIRCUMF> const& order) {
    uint cost = 0;
    for (uint pos=0; pos<CIRCUMF; pos++) {
        for (uint pos2=0; pos2<pos; pos2++) {
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
std::string BetaGrip::SimulatedAnnealing(
    uint nIter, double tempInit, double tempCool, ulong rseed
) {
    // initialise ordering
    auto order = std::array<uint, CIRCUMF>();
    auto orderBest = order;
    auto result = std::string(CIRCUMF, ' ');
    for (uint i=0; i<CIRCUMF; i++) {
        order[i] = i;
    }
    // Fisher-Yates shuffle for random initialisation
    rk_state rstate;
    rk_seed(rseed, &rstate);
    FYShuffle(order, CIRCUMF, rstate);
    auto cost = EvalCost(order);
    auto costBest = cost;

    double temp = tempInit;
    for (uint iter=0; iter<nIter; iter++) {
        auto orderNew = order;

        // swap a single pair each time
        {
            uint i = rk_interval(CIRCUMF-1, &rstate);
            uint j;
            do {
                j = rk_interval(CIRCUMF-1, &rstate);
            } while (j == i);
            uint temp = orderNew[i];
            orderNew[i] = orderNew[j];
            orderNew[j] = temp;
        }

        // accept new ordering if better
        // or if random acceptance function satisfied
        uint costNew = EvalCost(orderNew);
        bool update = false;
        if (costNew < cost) {
            update = true;
        } else {
            double costDelta = costNew - cost; // no worry about sign as delta>=0
            double rand = rk_double(&rstate);
            if (exp(-costDelta/temp) > rand) {
                update = true;
            }
        }
        // update ordering depending on above
        if (update) {
            order = orderNew;
            cost = costNew;
        }
        // save best found so far
        if (cost < costBest) {
            result = Idxs2String(order);
            costBest = cost;
        }
        // decay temperature
        temp *= tempCool;
    }
    return result;
}


std::string BetaGrip::GeneticEvolution(
    uint nGens, uint nPopu, uint nElite, uint nMerit, ulong rseed
) {
    auto result = std::string(CIRCUMF, ' ');
    uint costBest = -1;

    // initialise population
    auto geneBase = std::array<uint, CIRCUMF>();
    for (uint i=0; i<CIRCUMF; i++) {
        geneBase[i] = i;
    }
    rk_state rstate;
    rk_seed(rseed, &rstate);
    auto population = std::vector<std::array<uint, CIRCUMF>>(nPopu);
    auto ranked = std::vector<uint>(nPopu);
    for (uint i=0; i<nPopu; i++) {
        auto gene = geneBase;
        FYShuffle(gene, CIRCUMF, rstate);
        population[i] = std::move(gene);
        ranked[i] = i;
#if VERBOSE
        std::cout << "initial population:\n";
        auto str = Idxs2String(gene);
#endif
    }

    // do evolution
    for (uint gen=0; gen<nGens; gen++) {
        // find fitnesses and rank
        auto costs = std::vector<uint>(nPopu);
        for (uint i=0; i<nPopu; i++) {
            costs[i] = EvalCost(population[i]);
        }
        std::sort(ranked.begin(), ranked.end(), [&costs](uint i, uint j) {
            return costs[i] < costs[j];
        });
        // select top nElite for elitism
        auto survived = std::vector<bool>(nPopu, false);
        for (uint i=0; i<nElite; i++) {
            survived[ranked[i]] = true;
        }
        // select with tournament (binary)
        auto Revive = [&]() {
            uint player;
            do {
                player = rk_interval(nPopu, &rstate);
            } while (survived[player]);
            survived[player] = true;
            return player;
        };
        for (uint i=0; i<nMerit; i++) {
            auto p1 = Revive();
            auto p2 = Revive();
            auto loser = costs[p1]>costs[p2]? p1:p2;
            survived[loser] = false;
        }
        // breed using
        // OX1 from "Learning Bayesian Network Structures by searching
        // for the best ordering with genetic algorithms", 1996
        auto survivors = std::vector<uint>();
        survivors.reserve(nElite+nMerit);
        for (uint i=0; i<nPopu; i++) {
            if (survived[i]) {
                survivors.push_back(i);
                // std::cout << costs[i] << '\t';
                // auto str = Idxs2String(population[i]);
            }
        }
        break;
        // mutate
    }
    return result;
}