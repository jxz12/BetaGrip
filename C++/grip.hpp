#ifndef GRIP_HPP
#define GRIP_HPP

#include <string>
#include <vector>
#include <array>
#include <unordered_map>

#include "randomkit.h"

#define CIRCUMF 32
#define uint unsigned int
#define ulong unsigned long

/*
places letters around a wheel
e.g. 
           a b 
         h     c
         g     d
           f e

but in an optimal order that minimises
movement, defined as movements around
the wheel.
*/

class BetaGrip {
public:
    BetaGrip(const std::string& textPath);  // TODO: multiple text files
    std::string BruteForce();
    std::string SimulatedAnnealing(uint nIter, double tempInit, double tempCool, ulong rseed=0);
    std::string GeneticEvolution(uint nGens, uint nPopu, uint nElite, uint nMerit, ulong rseed=0);

private:
    std::unordered_map<char, uint> char2freq;
    std::unordered_map<char, uint> char2idx;
    std::array<char, CIRCUMF> idx2char;
    std::array<uint, CIRCUMF*CIRCUMF> freqMatrix;

    inline uint Two2OneD(uint const& row, uint const& col);
    inline std::string Idxs2String(std::array<uint, CIRCUMF> const& order);
    inline uint EvalCost(std::array<uint, CIRCUMF> const& order);
};

template<typename T>
void FYShuffle(T& arr, uint size, rk_state& rstate) {
    for (uint i=size-1; i>0; i--)
    {
        uint j = rk_interval(i, &rstate);
        auto temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
}

#endif