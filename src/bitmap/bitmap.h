#ifndef SNAPDRAGON_BITMAP_H
#define SNAPDRAGON_BITMAP_H
#include <vector>
#include "bitvector.h"
using namespace std;

/*
class BitmapIndex {
public:
    BitmapIndex();
    void loadIndex(const char* fname);
    void saveIndex(const char* fname);
protected:
    // every type of bitmap index has some BitVectors
    // total number of values indexed
    size_t                size;
    // the maximum number of values to hold in memory at once
    static const size_t   chunkSize=1000000;
    size_t                currentRow; // uncompressed position
    const char*           fname;
};
*/

// The BitSlicedIndex has one bitvector for each bit position
// The values may span more than one word but are expected
// to all be the same size.
// Compression is best when the values have been sorted.
template <class T>
class BitSlicedIndex {// : protected BitmapIndex {
public:
    // start a new bit sliced index
    BitSlicedIndex(int nwords, char* fname);
    ~BitSlicedIndex() {
        for(int i=0;i<nwords*nbits;i++)
            delete bvec[i];
    }
    void save() {}
    // append the referenced value
    void append(T *value);
    // reconstruct the value stored at position idx, return false if idx is invalid
    bool decode(size_t idx, T *value);
    
private:
    BitVector<T> **bvec;
    int nwords;
    char*           fname;
    T nbits;
    T **buffer; // square matrix of row oriented values (interface with user) or column oriented values (interface with bvec)
    size_t bufferCapacity; //= sizeof(T)*8; // might have to initialize in constructor
    size_t bufferOffset; // position of the beginning of the Buffer
    void transpose(uint64_t A[64]);
    void transpose(uint32_t A[32]);
    void transpose(uint16_t A[16]);
    // in constructor: allocate space for nwords vectors then for(i=0;i<nwords;i++) buffer[i].reserve(bufferCapacity);
    // when appending values just do rowBuffer[i].push_back(value[i]);
    // if the rowBuffer is full, transpose into uncompressed bitvectors into colBuffer and
    // invoke a bvec function to appendBits(T bits) on each value in the colBuffer
    // N.B. I dont care how bvec.appendBits() works, it might have it's own buffer of uncompressed words...
    // clear the vectors
    // The first time decode() is called, populate the colBuffer with sizeof(T)*8 uncompressed bits from each bvec
    // 
};
class RangeEncodedIndex {
public:
    RangeEncodedIndex(vector<uint32_t> &values, char *fname);
    ~RangeEncodedIndex() {
        for(int i=0;i<ranges.size();i++)
            delete bvec[i];
    }
    void save() {}
private:
    vector<uint32_t> ranges;
    BitVector<uint64_t> **bvec; // use 32bit words in the range encoded bitvector
    char *fname;
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
