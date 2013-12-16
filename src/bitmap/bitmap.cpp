// #include "bitmap.h"
// void BitmapIndex::loadIndex(const char* fname) {
//     
// }
// void BitmapIndex::saveIndex(const char* fname) {
//     
// }

template <class T>
BitSlicedIndex<T>::BitSlicedIndex(int nwords, const char* fname) {
    
    this->nwords = nwords;
    this->fname = fname;
    this->nbits = 8*sizeof(T);
    int nVectors = nwords*nbits;
    bvec = (BitVector<T>*) malloc(nVectors*sizeof(int*));
    for(int i=0;i<nVectors;i++)
        bvec[i] = new BitVector<T>();
    
    this->bufferCapacity = nbits;
    this->bufferOffset=0;

    this->buffer = (vector<T>*) malloc(nwords*sizeof(vector<T>));
    for(int i=0;i<nwords;i++)
        buffer[i].resize(bufferCapacity);
}

void transpose(unsigned long A[64]) {
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
void transpose(unsigned int A[32]) {
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
void transpose(unsigned short A[16]) {
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
            bvec[i].appendWord(buffer[i]);
        }
        bufferOffset = 0;
    }
}

template <class T>
bool BitSlicedIndex<T>::decode(size_t idx, T *value) {
    
}

// template <class T>
// RangeEncodedIndex<T>::RangeEncodedIndex(vector<T> vec, const char* fname) {
//     
// }
// 
// template <class T>
// void RangeEncodedIndex<T>::append(T value) {
//     
// }
// 
// template <class T>
// T RangeEncodedIndex<T>::decode(size_t idx) {
//     
// }

