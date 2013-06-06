#ifndef SNAPDRAGON_BITMAP_H
#define SNAPDRAGON_BITMAP_H
#include <vector>
#include "../bvec/bvec.h" // Bitvector class

class BitmapIndex {
public:
    void loadIndex(const char* fname);
    void saveIndex(const char* fname);
protected:
    // every type of bitmap index has some Bitvectors
    std::vector<Bitvector*> bvec;
    // each Bitvector represents a non-negative integer relative to some origin
    Bitvector*              bvecValues;
    int                     bvecOrigin;
    // total number of values indexed
    uint64_t                size;
private:
    // the maximum number of values to hold in memory at once
    static const uint32_t   chunkSize=1000000;
};

// The BitSlicedIndex has one bitvector for each bit position
// The values may span more than one word but are expected
// to all be the same size.
// Compression is best when the values have been sorted.
template <class T>
class BitSlicedIndex : public BitmapIndex {
public:
    // start a new bit sliced index
    BitSlicedIndex(unsigned int nwords, const char* fname);
    // append the given or referenced value
    void append(T value); // only works when nwords==1
    void append(T *value);
    // reconstruct the value stored at position idx
    T decode(size_t idx);
    T* decode(size_t idx);
    
private:
    uint32_t nwords;
};

// RangeEncodedIndex has one bitvector for each distinct
// value. The bitvector corresponding to x has a 1 in each
// position where the stored value >= x
template <class T>
class RangeEncodedIndex : public BitmapIndex {
public:
    // construct range encoding from a vector of numbers
    RangeEncodedIndex(vector<T> vec, const char* fname);
    // append the given value
    void append(T value);
    // reconstruct the value store at position idx
    T decode(size_t idx);
    
private:
};

#endif
