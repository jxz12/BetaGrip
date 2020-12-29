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

#include "peach.hpp"
#include "randomkit.h"

#define VERBOSE 1

inline uint PeachWheel::Two2OneD(uint const& row, uint const& col) {
    return row*CIRCUMF + col;
}

PeachWheel::PeachWheel(const std::string& textPath) {
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
    std::sort(idx2charAll.begin(), idx2charAll.end(), [&](char a, char b) {
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

inline std::string PeachWheel::Idxs2String(std::array<uint, CIRCUMF> const& order) {
    auto str = std::string(CIRCUMF, ' ');
    for (uint i=0; i<CIRCUMF; i++) {
        char c = idx2char[order[i]];
        str[i] = c;
    }
#if VERBOSE
    for (uint i=0; i<CIRCUMF; i++) {
        std::cout << (str[i]=='\n'? '_':str[i]);
    }
    std::cout << "\tc=" << EvalCost(order) << std::endl;
#endif
    return str;
}

// an exhaustive search of all permutations
std::string PeachWheel::BruteForce() {
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

uint PeachWheel::EvalCost(std::array<uint, CIRCUMF> const& order) {
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

std::string PeachWheel::RandomSample(uint nIter, ulong rseed) {
    auto order = std::array<uint, CIRCUMF>();
    for (uint i=0; i<CIRCUMF; i++) {
        order[i] = i;
    }
    auto result = std::string();
    uint best = -1;  // underflows to max
    rk_state rstate;
    rk_seed(rseed, &rstate);
    for (uint i=0; i<nIter; i++) {
        FYShuffle(order, CIRCUMF, rstate);
        auto cost = EvalCost(order);
        if (cost < best) {
            result = Idxs2String(order);
            best = cost;
        }
    }
    return result;
}

// simulated annealing
std::string PeachWheel::SimulatedAnnealing(
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
            std::swap(orderNew[i], orderNew[j]);
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


// genetic algorithm
// selection is performed with binary tournament
// breeding by ordered crossover
// mutation by chance of single swap per child
std::string PeachWheel::GeneticEvolution(
    uint nGens, uint nPopu, uint nElite, uint nMerit, double probMut, ulong rseed
) {
    rk_state rstate;
    rk_seed(rseed, &rstate);

    // initialise population
    auto geneBase = std::array<uint, CIRCUMF>();
    for (uint i=0; i<CIRCUMF; i++) {
        geneBase[i] = i;
    }
    auto population = std::vector<std::array<uint, CIRCUMF>>(nPopu, geneBase);
    auto ranked = std::vector<uint>(nPopu);
    for (uint i=0; i<nPopu; i++) {
        FYShuffle(population[i], CIRCUMF, rstate);
        ranked[i] = i;
    }
    // find fitnesses and rank them
    auto result = std::string(CIRCUMF, ' ');
    auto costs = std::vector<uint>(nPopu);
    uint costBest = -1;  // underflows to max
    auto Rank = [&]() {
        for (uint i=0; i<nPopu; i++) {
            costs[i] = EvalCost(population[i]);
            if (costs[i] < costBest) {
                result = Idxs2String(population[i]);
                costBest = costs[i];
            }
        }
        std::sort(ranked.begin(), ranked.end(), [&costs](uint i, uint j) {
            return costs[i] < costs[j];
        });
    };
    Rank();

    // do evolution
    for (uint gen=0; gen<nGens; gen++) {
        // select top nElite for elitism
        auto survivors = std::vector<bool>(nPopu, false);
        auto parents = std::vector<uint>();
        parents.reserve(nElite+nMerit);
        auto populationNew = std::vector<std::array<uint, CIRCUMF>>();
        populationNew.reserve(nPopu);
        for (uint i=0; i<nElite; i++) {
            uint elite = ranked[i];
            parents.push_back(elite);
            survivors[elite] = true;
            populationNew.push_back(population[elite]); // ensures best always kept
        }
        // select with tournament (binary)
        auto Sample = [&]()->uint {
            uint player;
            do {
                player = rk_interval(nPopu-1, &rstate);
            } while (survivors[player]);
            survivors[player] = true;
            return player;
        };
        for (uint i=0; i<nMerit; i++) {
            uint p1 = Sample();
            uint p2 = Sample();
            uint winner = costs[p1]<costs[p2]? p1:p2;
            uint loser = costs[p1]<costs[p2]? p2:p1;
            survivors[loser] = false;
            parents.push_back(winner);
        }
        // breed using
        // OX1 from "Learning Bayesian Network Structures by searching
        // for the best ordering with genetic algorithms", 1996
        for (uint i=nElite; i<nPopu; i++) {
            // select 2 parents
            uint mumIdx = parents[rk_interval(nMerit-1, &rstate)];
            uint dadIdx = parents[rk_interval(nMerit-1, &rstate)];
            auto mum = population[mumIdx];
            auto dad = population[dadIdx];
            // breed
            uint mumStart = rk_interval(CIRCUMF-1, &rstate);
            uint mumEnd = rk_interval(CIRCUMF-1, &rstate);
            if (mumEnd < mumStart) {
                std::swap(mumStart, mumEnd);
            }
            // copy selection from mum
            auto fromMum = std::vector<bool>(CIRCUMF, false);
            auto child = std::array<uint, CIRCUMF>();
            for (uint pos=mumStart; pos<=mumEnd; pos++) {
                auto inherited = mum[pos];
                child[pos] = inherited;
                fromMum[inherited] = true;
            }
            // fill in rest from dad
            uint dadPos = -1;
            for (uint pos=0; pos<CIRCUMF; pos++) {
                if (pos < mumStart || pos > mumEnd) {
                    while (fromMum[dad[++dadPos]]);  // lol
                    child[pos] = dad[dadPos];
                }
            }
            // mutate with probability
            double rand = rk_double(&rstate);
            if (rand <= probMut) {
                uint i = rk_interval(CIRCUMF-1, &rstate);
                uint j;
                do {
                    j = rk_interval(CIRCUMF-1, &rstate);
                } while (j == i);
                std::swap(child[i], child[j]);
            }
            populationNew.push_back(child);
        }
        population = std::move(populationNew);
        Rank();
    }
    return result;
}