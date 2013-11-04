#ifndef SNAPDRAGON_BITMAP_H
#define SNAPDRAGON_BITMAP_H
#include <vector>
#include "bitvector.h"
using namespace std;

// The BitSlicedIndex has one bitvector for each bit position
// The values may span more than one word but are expected
// to all be the same size.
// Compression is best when the values have been sorted.
template <class T>
class BitSlicedIndex {
public:
    // start a new index
    BitSlicedIndex(int nwords);
    // load an index from a file
    BitSlicedIndex(char* fname);
    // clean up
    ~BitSlicedIndex() {
        for(int i=0;i<nwords*nbits;i++)
            delete bvec[i];
    }
    // save an index to the file
    void saveIndex(char* fname);
    // append the referenced value
    void append(T *value);
    void append(T *value, int x);
    // reconstruct the value stored at position idx, return false if idx is invalid
    bool decode(size_t idx, T *value);
    size_t size() {return nValues;}
    

private:
    BitVector<T> **bvec;
    int nwords;
    T nbits; // 8*sizeof(T)
    size_t nValues; // total number of values represented - same as size of each bvec
    T **buffer; // square matrix of row oriented values (interface with user) or column oriented values (interface with bvec)
    size_t bufferCapacity; //= sizeof(T)*8; had to initialize in constructor
    size_t bufferOffset; // offset within the buffer - used in append()
    size_t bufferStart; // position of the beginning of the Buffer
    void fillBuffer(size_t idx); // fill the buffer with words to cover idx, returns false if idx too large

    void transpose(uint64_t A[64]);
    void transpose(uint32_t A[32]);
    void transpose(uint16_t A[16]);
};
class RangeEncodedIndex {
public:
    RangeEncodedIndex(vector<uint32_t> &values);
    RangeEncodedIndex(char *fname);
    ~RangeEncodedIndex() {
        for(int i=0;i<ranges.size();i++)
            delete bvec[i];
    }
    void saveIndex(char *fname);
    // reconstruct the value stored at position idx, return false if idx is invalid
    bool decode(size_t idx, uint32_t *value);
    size_t size() {return nValues;}
private:
    vector<uint32_t> ranges;
    BitVector<uint64_t> **bvec;
    size_t nValues;
    size_t bufferStart;
    size_t bufferSize;
    uint32_t *buffer, *buffer0;
    void fillBuffer(size_t idx);
};


#include "bitmap.tpp"
/*
// RangeEncodedIndex has one bitvector for each distinct
// value. The bitvector corresponding to x has a 1 in each
// position where the stored value >= x
template <class T>
class RangeEncodedIndex : protected BitmapIndex {
public:
    // construct range encoding from a vector of numbers
    RangeEncodedIndex(vector<T> vec, const char* fname);
    // append the given value
    void append(T value);
    // reconstruct the value store at position idx
    T decode(size_t idx);
    
private:
    vector<BitVector<T>*> bvec;
};
*/


#endif
