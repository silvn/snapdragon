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
    void inflateWord(T *word, T wordStart); // fills a word sized uncompressed bitvector starting at wordStart
    void appendFill0(T length); // extend the bitvector by length 0 bits
    void appendFill1(T length); // extend the bitvector by length 1 bits
    bool operator [](T pos); // check the bit at this position
    T nextOne(T pos); // returns the position of the next set bit
    void fillSetBits(T wordStart, uint32_t *A, uint32_t v);

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
    T activeWordEnd;   // uncompressed bit position of last bit in a word
    void seek(T wordStart);
};
inline int ffs(unsigned long long bits) { return __builtin_ffsll(bits); }
inline int ffs(unsigned long bits)      { return __builtin_ffsl (bits); }
inline int ffs(unsigned int bits)       { return __builtin_ffs  (bits); }
#include "bitvector.tpp"
#endif
