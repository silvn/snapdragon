#ifndef SNAPDRAGON_LCBitSlicedIndex_H
#define SNAPDRAGON_LCBitSlicedIndex_H
#include <vector>
#include "bitvector.hpp"
using namespace std;

// The LCBitSlicedIndex is similar to the LCBitSlicedIndex class
// where the indexed values have relatively low cardinality
// BitSlicedIndex is used for high cardinality data and lacks the
// vectors of distinct values (dValues) and associated counts (dFrequency).
template <class T>
class LCBitSlicedIndex {
public:
    // start a new index
    LCBitSlicedIndex();
    // load an index from a file
    LCBitSlicedIndex(char* fname);
    // clean up
    ~LCBitSlicedIndex() {
        for(int i=0;i<nbits;i++)
            delete bvec[i];
    }
    // save an index to the file
    void saveIndex(char* fname);
    // append the referenced value
    void append(T *value);
    // reconstruct the value stored at position idx, return false if idx is invalid
    bool decode(size_t idx, T *value);
    size_t size() {return nValues;}

    // query functions - supply a minimal min or a maximal max for one sided range query
    BitVector<T>* continuousRange(const T min, const T max);
    // given a list of values, convert into a set of continuousRange queries
    // and OR the result vectors
    BitVector<T>* discreteRange(vector<T> vals);        

private:
    BitVector<T> **bvec;
    T nbits; // 8*sizeof(T)
    size_t nValues; // total number of values represented - same as size of each bvec
    vector<T> dValues; // vector of distinct sorted values
    vector<uint32_t> dFrequency; // corresponding vector of frequencies
    
    uint64_t bv[256]; // 256*64 bits = 2^14 = 16384
    uint32_t counts[16384];
    vector<T> overflow;
    void insertValue(T x); // called by append
    void makeDistinct();   // called by saveIndex

    BitVector<T>* evalGE(T x);

    T *buffer; // square matrix of row oriented values (interface with user) or column oriented values (interface with bvec)
    size_t bufferCapacity; //= sizeof(T)*8; had to initialize in constructor
    size_t bufferOffset; // offset within the buffer - used in append()
    size_t bufferStart; // position of the beginning of the Buffer
    void fillBuffer(size_t idx); // fill the buffer with words to cover idx, returns false if idx too large

    void transpose(uint64_t A[64]);
    void transpose(uint32_t A[32]);
    void transpose(uint16_t A[16]);
};


template <class T>
LCBitSlicedIndex<T>::LCBitSlicedIndex() {
    
    this->nbits = 8*sizeof(T);
    int nVectors = nbits;
    bvec = (BitVector<T>**) malloc(nVectors*sizeof(BitVector<T>*));
    for(int i=0;i<nVectors;i++)
        bvec[i] = new BitVector<T>();
    nValues=0;

    memset(bv,0,256*sizeof(uint64_t));
    memset(counts,0,16384*sizeof(uint32_t));
    
    this->bufferCapacity = nbits;
    this->bufferOffset=0;
    this->bufferStart=0;

    this->buffer = (T*) malloc(bufferCapacity*sizeof(T*));
}


/*
    // write the number of values (length)
    fwrite(&nValues, sizeof(size_t),1,fp);
    // write the number of distinct values
    size_t nDistinct = dValues.size();
    fwrite(&nDistinct, sizeof(size_t),1,fp);
    // write the distinct values
    fwrite(dValues.data(), sizeof(T), nDistinct, fp);
    // write the associated frequencies
    fwrite(dFrequency.data(), sizeof(uint32_t), nDistinct, fp);
    
*/
template <class T>
LCBitSlicedIndex<T>::LCBitSlicedIndex(char* fname) {
    this->nbits = 8*sizeof(T);
    FILE *fp;
    fp = fopen(fname, "rb");
    // read number of values (length)
    size_t result = fread(&nValues, sizeof(size_t),1,fp);
    if (result != 1) {fputs ("Reading error2",stderr); exit (3);}
    // read number of distinct values
    size_t nDistinct;
    result = fread(&nDistinct, sizeof(size_t),1,fp);
    if (result != 1) {fputs ("Reading error2",stderr); exit (3);}
    // read distinct values and corresponding frequencies
    dValues.resize(nDistinct);
    result = fread(dValues.data(), sizeof(T), nDistinct, fp);
    if (result != nDistinct) {fputs ("Reading error2",stderr); exit (3);}
    dFrequency.resize(nDistinct);
    result = fread(dFrequency.data(), sizeof(T), nDistinct, fp);
    if (result != nDistinct) {fputs ("Reading error2",stderr); exit (3);}
    // read each bitvector (make this optional for when using histogram functions)
    int nVectors = nbits;
    bvec = (BitVector<T>**) malloc(nVectors*sizeof(BitVector<T>*));
    for(int i=0;i<nVectors;i++) {
        T nT,*bvbuffer;
        result = fread(&nT,sizeof(int),1,fp);
        if (result != 1) {fputs ("Reading error3",stderr); exit (3);}
        bvbuffer = (T*) malloc(nT*sizeof(T));
        result = fread(bvbuffer,sizeof(T), nT, fp);
        bvec[i] = new BitVector<T>(bvbuffer);
        free(bvbuffer);
    }
    fclose(fp);
     
    this->bufferCapacity = nbits;
    this->bufferOffset=0;
    this->bufferStart=0;

    this->buffer = (T*) malloc(bufferCapacity*sizeof(T));
    fillBuffer(0);
}

template <class T>
void LCBitSlicedIndex<T>::transpose(uint64_t A[64]) {
    int j, k;
    uint64_t m, t;
    m = 0x00000000FFFFFFFFULL;
    for (j = 32; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 64; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}
template <class T>
void LCBitSlicedIndex<T>::transpose(uint32_t A[32]) {
    int j, k;
    uint32_t m, t;
    
    m = 0x0000FFFF;
    for (j = 16; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 32; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}
template <class T>
void LCBitSlicedIndex<T>::transpose(uint16_t A[16]) {
    int j, k;
    uint16_t m, t;
    
    m = 0x00FF;
    for (j = 8; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 16; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}

template <class T>
void LCBitSlicedIndex<T>::insertValue(T x) {
    if (x < 16364) {
        bv[x >> 6] |= (1ULL << (x & 63));
        counts[x]++;
    }
    else
        overflow.push_back(x);
}

template <class T>
void LCBitSlicedIndex<T>::makeDistinct() {
    T v;
    for(int i=0;i<256;i++) {
        uint64_t bits = bv[i];
        while (bits) {
            v = i*64 + __builtin_ffsll(bits) - 1;
            dValues.push_back(v);
            dFrequency.push_back(counts[v]);
            bits &= bits-1;
        }
    }
    if (overflow.size() > 0) {
        sort(overflow.begin(),overflow.end());
        // iterate over the sorted values
        dValues.push_back(overflow[0]);
        int last=dFrequency.size();
        dFrequency.push_back(1);
        for(int i=1;i<overflow.size();i++) {
            if (overflow[i] == dValues.back())
                dFrequency[last]++;
            else {
                dValues.push_back(overflow[i]);
                dFrequency.push_back(1);
                last++;
            }
        }
    }
}


template <class T>
void LCBitSlicedIndex<T>::append(T* value) {
    insertValue(*value);
    buffer[bufferOffset] = *value;
    bufferOffset++;
    if (bufferOffset == nbits) {
        this->transpose(buffer);
        for(int j=0;j<nbits;j++)
            bvec[j]->appendWord(buffer[j]);
        bufferOffset = 0;
    }
    nValues++;
}


template <class T>
void LCBitSlicedIndex<T>::fillBuffer(size_t idx) {
    // fprintf(stderr,"fillBuffer(%zi)\n",idx);
    // in each bitvector, fetch an uncompressed bitvector word that contains idx
    // populate the buffers by transposing the uncompressed words
    bufferStart = idx & ~(nbits - 1);
    for(int j=0;j<nbits;j++)
        bvec[j]->inflateWord(buffer+j,bufferStart);
    this->transpose(buffer);
}

template <class T>
bool LCBitSlicedIndex<T>::decode(size_t idx, T *value) {
    // fprintf(stderr,"decode(%zi) nValues: %zi bufferStart: %zi nbits: %zi\n",idx,nValues,bufferStart,nbits);
    if (idx >= nValues) return false;
    if ((idx >= bufferStart + nbits) || idx < bufferStart)
        fillBuffer(idx);
    
    int offset= idx - bufferStart;
    *value = buffer[offset];
    return true;
}

template <class T>
void LCBitSlicedIndex<T>::saveIndex(char *fname) {
    // flush the buffer
    if (bufferOffset > 0) {
        // fill buffer with 0's
        for(int j=bufferOffset;j<nbits;j++)
            buffer[j] = (T)0;
        this->transpose(buffer);
        for(int j=0;j<nbits;j++)
            bvec[j]->appendWord(buffer[j]);
    }
    // populate the dValues vector
    makeDistinct();
    // open an output file
    FILE *fp;
    fp = fopen(fname, "wb");
    // write the number of values (length)
    fwrite(&nValues, sizeof(size_t),1,fp);
    // write the number of distinct values
    size_t nDistinct = dValues.size();
    fwrite(&nDistinct, sizeof(size_t),1,fp);
    // write the distinct values
    fwrite(dValues.data(), sizeof(T), nDistinct, fp);
    // write the associated frequencies
    fwrite(dFrequency.data(), sizeof(uint32_t), nDistinct, fp);
    // write each bitvector
    for(int i=0; i<nbits; i++) {
        T *buf;
        size_t nT = bvec[i]->dump(&buf);
        fwrite(&nT,sizeof(int),1,fp);
        fwrite(buf,sizeof(T),nT,fp);
        free(buf);
    }
    fclose(fp);
}

inline int ffs(unsigned long long bits) { return __builtin_ffsll(bits); }
inline int ffs(unsigned long bits)      { return __builtin_ffsl (bits); }
inline int ffs(unsigned int bits)       { return __builtin_ffs  (bits); }
inline int clz(unsigned long long bits) { return __builtin_clzll(bits); }
inline int clz(unsigned long bits)      { return __builtin_clzl (bits); }
inline int clz(unsigned int bits)       { return __builtin_clz  (bits); }

// create a bitvector marking rows where the value >= x
template <class T>
BitVector<T>* LCBitSlicedIndex<T>::evalGE(T x) {
    if (x == dValues[0]) {
        // return a bitvector of all ones
        BitVector<T> *res = new BitVector<T>();
        res->appendFill1(nValues);
        return res;
    }

    BitVector<T>* results = new BitVector<T>();
    BitVector<T>* partial = new BitVector<T>();
    results->appendFill0(nValues);
    partial->appendFill1(nValues);
    int msb = clz(dValues.back());
    int lsb = nbits - ffs(x);
    T selector = (T)1 << (nbits-1);
    for (int b = msb; b <= lsb; b++)
        if (x & (selector >> b))
            partial &= bvec[b];
        else
            results |= partial & bvec[b]; 
    results |= partial;
    return results;
}

// query functions - supply a minimal min or a maximal max for one sided range query
template <class T>
BitVector<T>* LCBitSlicedIndex<T>::continuousRange(const T min, const T max) {
    // find lower and upper bound based on min and max
    typename vector<T>::iterator low,up;
    low=lower_bound (dValues.begin(), dValues.end(), min);
    if (low == dValues.end()) {
        // return empty bitvector
        BitVector<T> *res = new BitVector<T>();
        res->appendFill0(nValues);
        return res;
    }
    up= upper_bound (dValues.begin(), dValues.end(), max);

    // >= *low _AND_ NOT >= *up
    BitVector<T> *gelow = evalGE(*low);
    BitVector<T> *geup  = evalGE(*up);
    geup->flip();
    return gelow & geup;
}

// given a list of values, convert into a set of continuousRange queries
// and OR the result vectors
template <class T>
BitVector<T>* LCBitSlicedIndex<T>::discreteRange(vector<T> vals) {
    vector<int> ranges;
    typename vector<T>::iterator low;
    for(typename vector<T>::iterator it=vals.begin(); it != vals.end(); it++) {
        low = lower_bound(dValues.begin(), dValues.end(), *it);
        if (*low == *it)
            ranges.push_back(distance(low,dValues.begin()));
    }
    BitVector<T>* result = new BitVector<T>();
    result->appendFill0(nValues);
    if (ranges.size() == 0)
        return result;
    // ranges is a set of offsets into dValues
    sort(ranges.begin(),ranges.end());
    int from = ranges[0];
    int to = from;
    for(int i=1;i<ranges.size();i++) {
        if (ranges[i] > to+1) {
            result |= continuousRange(dValues[from], dValues[to]);
            from = ranges[i];
        }
        to = ranges[i];
    }
    result |= continuousRange(dValues[from], dValues[to]);
    return result;
}



#endif