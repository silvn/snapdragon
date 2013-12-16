#ifndef SNAPDRAGON_RANGEENCODEDINDEX_H
#define SNAPDRAGON_RANGEENCODEDINDEX_H
#include <vector>
#include "bitvector.hpp"
using namespace std;

// RangeEncodedIndex has one bitvector for each distinct
// value. The bitvector corresponding to x has a 1 in each
// position where the stored value >= x
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



RangeEncodedIndex::RangeEncodedIndex(vector<uint32_t> &values) {
    nValues = values.size();
    vector<uint32_t>::iterator it;
    // find the distinct values and store in ranges vector
    if (0) {
        ranges = values;
        sort(ranges.begin(),ranges.end());
        it = unique(ranges.begin(),ranges.end());
        ranges.resize(distance(ranges.begin(),it));
    }
    if (1) {
        ranges.clear();
        // iterate over the values once, using a small uncompressed bitvector to mark the values < 256
        // push larger values into an overflow vector to sortuniq later
        uint64_t bv[256]; // 256*64 bits = 2^14 = 16384
        memset(bv,0,256*sizeof(uint64_t));
        vector<uint32_t> overflow;
        for(it = values.begin(); it != values.end();it++) {
            if (*it < 16364) {
                // set the bit
                bv[*it >> 6] |= (1ULL << (*it & 63));
            }
            else {
                overflow.push_back(*it);
            }
        }
        for(int i=0;i<256;i++) {
            uint64_t bits = bv[i];
            while (bits) {
                ranges.push_back(i*64 + __builtin_ffsll(bits) - 1);
                bits &= bits-1;
            }
        }
        if (overflow.size() > 0) {
            sort(overflow.begin(),overflow.end());
            it = unique(overflow.begin(), overflow.end());
            overflow.resize(distance(overflow.begin(),it));
            ranges.insert(ranges.end(), overflow.begin(), overflow.end());
        }
    }
    
    // populate a bvec for each distinct value
    bvec = (BitVector<uint64_t>**) malloc(ranges.size()*sizeof(BitVector<uint64_t>*));
    for(int i=0;i<ranges.size();i++)
        bvec[i] = new BitVector<uint64_t>();

    // iterate over the values again and append to each of the relevant bitvectors
    // compare the bitvector index of the current value to the previous value and 
    // for each index between them (where a run of 1's or 0's ends), appendFill[01](runLength)
    
    // we need an array to hold the start positions of each run
    uint32_t *runStart;
    runStart = (uint32_t*) calloc(ranges.size(),sizeof(uint32_t)); // init to 0

    uint32_t currentPos = 0;
    it = values.begin();
    vector<uint32_t>::iterator lb = lower_bound(ranges.begin(), ranges.end(), *it);
    int prevIdx = lb - ranges.begin();
    it++;
    currentPos++;
    while(it != values.end()) {
        int idx;
        if (*it < ranges[prevIdx]) {
            lb = lower_bound(ranges.begin(), lb, *it);
            idx = lb - ranges.begin();
            for(int i=idx+1;i<=prevIdx;i++) {
                bvec[i]->appendFill1(currentPos - runStart[i]);
                runStart[i] = currentPos;
            }
        }
        else if (*it > ranges[prevIdx]) {
            lb = lower_bound(lb+1, ranges.end(), *it);
            idx = lb - ranges.begin();
            for(int i=prevIdx+1; i<=idx; i++) {
                bvec[i]->appendFill0(currentPos - runStart[i]);
                runStart[i] = currentPos;
            }
        }
        prevIdx = idx;
        it++;
        currentPos++;
    }
    // finish each bitvector
    for(int i=0;i<=prevIdx;i++)
        bvec[i]->appendFill1(currentPos - runStart[i]);
    for(int i=prevIdx+1;i<ranges.size();i++)
        bvec[i]->appendFill0(currentPos - runStart[i]);
    delete runStart;
}

RangeEncodedIndex::RangeEncodedIndex(char *fname) {
    FILE *fp;
    fp = fopen(fname, "rb");

    size_t nVectors;
    size_t result = fread(&nVectors,sizeof(size_t),1,fp);
    if (result != 1) {fputs ("Reading error",stderr); exit (3);}

    result = fread(&nValues,sizeof(size_t),1,fp);
    if (result != 1) {fputs ("Reading error",stderr); exit (3);}

    ranges.resize(nVectors);
    result = fread(ranges.data(), sizeof(uint32_t),nVectors,fp);
    if (result != nVectors) {fputs ("Reading error",stderr); exit (3);}

    bvec = (BitVector<uint64_t>**) malloc(nVectors*sizeof(BitVector<uint64_t>*));
    for(int i=0;i<nVectors;i++) {
        uint64_t *bvbuffer;
        size_t nWords;
        result = fread(&nWords,sizeof(size_t),1,fp);
        if (result != 1) {fputs ("Reading error",stderr); exit (3);}
        bvbuffer = (uint64_t*) malloc(nWords*sizeof(uint64_t));
        result = fread(bvbuffer,sizeof(uint64_t), nWords, fp);
        bvec[i] = new BitVector<uint64_t>(bvbuffer);
        free(bvbuffer);
    }
    fclose(fp);
    bufferStart=0;
    bufferSize = 1 << 6; // tweak me
    buffer = (uint32_t*) malloc(bufferSize * sizeof(uint32_t));
    buffer0 = (uint32_t*) malloc(bufferSize * sizeof(uint32_t));
    for(int i=0;i<bufferSize;i++) buffer0[i] = ranges[0];
    fillBuffer(0);
}

void RangeEncodedIndex::saveIndex(char *fname) {
    // open an output file
    FILE *fp;
    fp = fopen(fname, "wb");
    // write the number of bitvectors
    size_t nVectors = ranges.size();
    fwrite(&nVectors,sizeof(size_t),1,fp);
    // write the number of values
    fwrite(&nValues,sizeof(size_t),1,fp);
    // write the distinct values
    fwrite(ranges.data(),sizeof(uint32_t),nVectors,fp);
    // write each bitvector
    for(int i=0;i<nVectors;i++) {
        uint64_t *buf;
        size_t nWords = bvec[i]->dump(&buf);
        fwrite(&nWords,sizeof(size_t),1,fp);
        fwrite(buf,sizeof(uint64_t),nWords,fp);
        free(buf);
        fprintf(stderr,"saveIndex(%s) %i/%zi nWords %zi\n",fname,i,nVectors,nWords);
    }
    fclose(fp);
}

void RangeEncodedIndex::fillBuffer(size_t idx) {
    bufferStart = idx & ~(bufferSize - 1);
    memcpy(buffer, buffer0, bufferSize*sizeof(uint32_t));
    uint64_t bits;
    for(int i=1;i<ranges.size();i++) {
        bvec[i]->inflateWord(&bits,bufferStart);
        if (bits == 0) break;
        while (bits) {
            buffer[__builtin_ffsll(bits)-1] = ranges[i];
            bits &= bits - 1ULL;
        }
    }
}

bool RangeEncodedIndex::decode(size_t idx, uint32_t *value) {
    if (idx >= nValues) return false;
    if ((idx >= bufferStart + bufferSize) || idx < bufferStart)
        fillBuffer(idx);
    
    int offset= idx - bufferStart;
    *value = buffer[offset];
    return true;
}
#endif