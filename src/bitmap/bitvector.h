#ifndef SNAPDRAGON_BITVECTOR_H
#define SNAPDRAGON_BITVECTOR_H
#include <vector>
using namespace std;
#define ZEROFILL 0
#define ONEFILL 1
#define LITERAL 2
#define NONEOFTHEABOVE 3

template <class T>
class BitVector {
public:
    BitVector(); // constructor
    BitVector(T *buf);
    size_t dump(T **buf);

    void appendWord(T word); // appends word shaped uncompressed bitvector
    bool nextWord(T *word, T *bit_pos); // fills a word shaped uncompressed bitvector, false if end of vector
    void appendFill0(T length); // extend the bitvector by length 0 bits
    void appendFill1(T length); // extend the bitvector by length 1 bits
    bool operator [](T pos); // check the bit at this position
    T nextOne(T pos); // returns the position of the next set bit
private:
    vector<T> words; // a mix of literal and fill words. MSB of fill words indicate the type of fill
    vector<T> isFill; // not using the vector<bool> specialization because we need to dump the data

    T count;   // cache the number of set bits
    T size;    // bits in the uncompressed bitvector
    T nbits;   // number of bits per LITERAL word
    T shiftby; // log base 2 of nbits

    // active word used for iterating
    
    T activeWordIdx;
    char activeWordType; // {ONEFILL, ZEROFILL, LITERAL}
    T activeWordStart; // uncompressed bit position at start of word
};
#include "bitvector.tpp"
#endif
