// #include "bitmap.h"
// void BitmapIndex::loadIndex(const char* fname) {
//     
// }
// void BitmapIndex::saveIndex(const char* fname) {
//     
// }
#include <bitset>

template <class T>
BitSlicedIndex<T>::BitSlicedIndex(int nwords, char* fname) {
    
    this->nwords = nwords;
    this->fname = fname;
    this->nbits = 8*sizeof(T);
    int nVectors = nwords*nbits;
    bvec = (BitVector<T>**) malloc(nVectors*sizeof(BitVector<T>*));
    for(int i=0;i<nVectors;i++)
        bvec[i] = new BitVector<T>();
    
    this->bufferCapacity = nbits;
    this->bufferOffset=0;

    this->buffer = (T**) malloc(nwords*sizeof(T*));
    for(int i=0;i<nwords;i++)
        buffer[i] = (T*) malloc(bufferCapacity*sizeof(T));
}

template <class T>
void BitSlicedIndex<T>::transpose(uint64_t A[64]) {
    int j, k;
    unsigned long m, t;
    
    m = 0x00000000FFFFFFFF;
    for (j = 32; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 64; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}
template <class T>
void BitSlicedIndex<T>::transpose(uint32_t A[32]) {
    int j, k;
    unsigned int m, t;
    
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
void BitSlicedIndex<T>::transpose(uint16_t A[16]) {
    int j, k;
    unsigned short m, t;
    
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
void BitSlicedIndex<T>::append(T* value) {
    for(int i=0; i<nwords; i++)
        buffer[i][bufferOffset] = value[i];
    bufferOffset++;
    if (bufferOffset == nbits) {
        for(int i=0;i<nwords;i++) {
            transpose(buffer[i]);
            for(int j=0;j<nbits;j++)
                bvec[i*nbits + j]->appendWord(buffer[i][j]);
        }
        bufferOffset = 0;
    }
}

template <class T>
bool BitSlicedIndex<T>::decode(size_t idx, T *value) {
    
}

RangeEncodedIndex::RangeEncodedIndex(vector<uint32_t> &values, char *fname) {
    this->fname = fname;
    // fprintf(stderr,"RangeEncodedIndex(%zi values, %s)\n",values.size(),fname);
    vector<uint32_t>::iterator it;
    // find the distinct values and store in ranges vector
    if(0) {
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
        for(it = values.begin(); it < values.end();it++) {
            if (*it < 16364) {
                // set the bit
                bv[*it >> 6] |= (1ULL << (*it & 63));
            }
            else {
                overflow.push_back(*it);
            }
        }
        sort(overflow.begin(),overflow.end());
        it = unique(overflow.begin(), overflow.end());
        overflow.resize(distance(overflow.begin(),it));
        for(int i=0;i<256;i++) {
            uint64_t bits = bv[i];
            while (bits) {
                ranges.push_back(i*64 + __builtin_ffsll(bits) - 1);
                bits &= bits-1;
            }
        }
        ranges.insert(ranges.end(), overflow.begin(), overflow.end());
    }
    // fprintf(stderr,"ranges.size() = %zi\n",ranges.size());
    // populate a bvec for each distinct value
    bvec = (BitVector<uint64_t>**) malloc(ranges.size()*sizeof(BitVector<uint64_t>*));
    for(int i=0;i<ranges.size();i++)
        bvec[i] = new BitVector<uint64_t>();

    // fprintf(stderr,"finished allocating bvec array\n");
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
    while(it < values.end()) {
        // fprintf(stderr,"ranges[%i]:%u values[%u]:%u values.size(): %zi\n",prevIdx,ranges[prevIdx],currentPos,*it, values.size());
        if (*it != ranges[prevIdx]) {
            lb = lower_bound(ranges.begin(), ranges.end(), *it);
            int idx = lb - ranges.begin();
            // fprintf(stderr,"idx:%i\n",idx);
            if (idx < prevIdx) {
                for(int i=idx+1;i<=prevIdx;i++) {
                    // fprintf(stderr,"bvec[%i]->appendFill1(%u)\n",i,currentPos-runStart[i]);
                    bvec[i]->appendFill1(currentPos - runStart[i]);
                    runStart[i] = currentPos;
                }
            }
            else {
                for(int i=prevIdx+1; i<=idx; i++) {
                    // fprintf(stderr,"bvec[%i]->appendFill0(%u)\n",i,currentPos-runStart[i]);
                    bvec[i]->appendFill0(currentPos - runStart[i]);
                    runStart[i] = currentPos;
                }
            }
            prevIdx = idx;
        }
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