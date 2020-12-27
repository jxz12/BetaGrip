#include <string>
#include <vector>
#include <array>
#include <unordered_map>

#define CIRCUMF 16

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

class GodOfKeyboards {
public:
    GodOfKeyboards(const std::string& textPath);  // TODO: multiple text files
    std::array<char, CIRCUMF> BranchAndBound();
    std::array<char, CIRCUMF> SimulatedAnnealing(int nIterations, float cooling);
    std::array<char, CIRCUMF> GeneticEvolution(int nIterations, float mutationProb);

private:
    std::unordered_map<char, unsigned> char2freq;
    std::unordered_map<char, unsigned> char2idx;
    std::array<char, CIRCUMF> idx2char;
    std::array<unsigned, CIRCUMF*CIRCUMF> freqMatrix;
    inline unsigned Two2OneD(unsigned row, unsigned col);
};