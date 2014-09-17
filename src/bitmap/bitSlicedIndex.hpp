#ifndef SNAPDRAGON_BITSLICEDINDEX_H
#define SNAPDRAGON_BITSLICEDINDEX_H
#include <vector>
#include "bitvector.hpp"
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
    // reconstruct the value stored at position idx, return false if idx is invalid
    bool decode(size_t idx, T *value);
    size_t size() {return nValues;}
    
    // search index for referenced value, if found set pos and return true
    bool find(T *value, size_t *pos);
    

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

    // void transpose(uint64_t A[64]);
    // void transpose(uint32_t A[32]);
    // void transpose(uint16_t A[16]);
};

template <class T>
BitSlicedIndex<T>::BitSlicedIndex(int nwords) {
    
    this->nwords = nwords;
    this->nbits = 8*sizeof(T);
    int nVectors = nwords*nbits;
    bvec = (BitVector<T>**) malloc(nVectors*sizeof(BitVector<T>*));
    for(int i=0;i<nVectors;i++)
        bvec[i] = new BitVector<T>();
    nValues=0;
    
    this->bufferCapacity = nbits;
    this->bufferOffset=0;
    this->bufferStart=0;

    this->buffer = (T**) malloc(nwords*sizeof(T*));
    for(int i=0;i<nwords;i++)
        buffer[i] = (T*) malloc(bufferCapacity*sizeof(T));
}

template <class T>
BitSlicedIndex<T>::BitSlicedIndex(char* fname) {
    this->nbits = 8*sizeof(T);
    // fprintf(stderr,"BitSlicedIndex(%s) nbits %zi\n",fname,nbits);
    // fread() the number of words
    FILE *fp;
    fp = fopen(fname, "rb");
    size_t result = fread(&nwords,sizeof(int),1,fp);
    if (result != 1) {fputs ("Reading error1",stderr); exit (3);}
    result = fread(&nValues, sizeof(size_t),1,fp);
    if (result != 1) {fputs ("Reading error2",stderr); exit (3);}
    int nVectors = nwords*nbits;
    bvec = (BitVector<T>**) malloc(nVectors*sizeof(BitVector<T>*));
    for(int i=0;i<nVectors;i++) {
        size_t nT;
        T *bvbuffer;
        result = fread(&nT,sizeof(size_t),1,fp);
        if (result != 1) {fputs ("Reading error3",stderr); exit (3);}
        bvbuffer = (T*) malloc(nT*sizeof(T));
        result = fread(bvbuffer,sizeof(T), nT, fp);
        bvec[i] = new BitVector<T>(bvbuffer);
        free(bvbuffer);
    }
    fclose(fp);
    
    // fprintf(stderr,"BitSlicedIndex(%s) nValues %zi nVectors %i\n",fname,nValues,nVectors);
 
    this->bufferCapacity = nbits;
    this->bufferOffset=0;
    this->bufferStart=0;

    this->buffer = (T**) malloc(nwords*sizeof(T*));
    for(int i=0;i<nwords;i++) {
        buffer[i] = (T*) malloc(bufferCapacity*sizeof(T));
        for(int j=0;j<nbits;j++)
            bvec[i*nbits+j]->inflateNextWord(buffer[i]+j,bufferStart);
        transpose(buffer[i]);
    }
//    fillBuffer(0);
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
    nValues++;
}


template <class T>
void BitSlicedIndex<T>::fillBuffer(size_t idx) {
    // in each bitvector, fetch an uncompressed bitvector word that contains idx
    // populate the buffers by transposing the uncompressed words

    // if this is the next buffer use a simple bitvector function
    if (idx == bufferStart + nbits) {
        bufferStart = idx;
        for(int i=0;i<nwords;i++) {
            for(int j=0;j<nbits;j++)
                bvec[i*nbits+j]->inflateNextWord(buffer[i]+j,bufferStart);
            transpose(buffer[i]);
        }
    }
    else { // assume random access
        // fprintf(stderr,"fillBuffer(%zi)\n",idx);
        bufferStart = idx & ~(nbits - 1);
        for(int i=0;i<nwords;i++) {
            for(int j=0;j<nbits;j++)
                bvec[i*nbits+j]->inflateWord(buffer[i]+j,bufferStart);
            transpose(buffer[i]);
        }
    }
}

template <class T>
bool BitSlicedIndex<T>::decode(size_t idx, T *value) {
    // fprintf(stderr,"decode(%zi) nValues: %zi bufferStart: %zi nbits: %zi\n",idx,nValues,bufferStart,nbits);
    if (idx >= nValues) return false;
    if ((idx >= bufferStart + nbits) || idx < bufferStart)
        fillBuffer(idx);
    
    int offset= idx - bufferStart;
    for(int i=0;i<nwords;i++)
        value[i] = buffer[i][offset];
    return true;
}

template <class T>
void BitSlicedIndex<T>::saveIndex(char *fname) {
    // flush the buffer
    if (bufferOffset > 0) {
        for(int i=0;i<nwords;i++) {
            // fill buffer with 0's
            for(int j=bufferOffset;j<nbits;j++)
                buffer[i][j] = (T)0;
            transpose(buffer[i]);
            for(int j=0;j<nbits;j++)
                bvec[i*nbits + j]->appendWord(buffer[i][j]);
        }
    }
    // open an output file
    FILE *fp;
    fp = fopen(fname, "wb");
    // write the number of words (width)
    fwrite(&nwords, sizeof(int),1,fp);
    // write the number of values (length)
    fwrite(&nValues, sizeof(size_t),1,fp);
    // write each bitvector
    for(int i=0; i<nwords*nbits; i++) {
        T *buf;
        size_t nT = bvec[i]->dump(&buf);
        fwrite(&nT,sizeof(size_t),1,fp);
        fwrite(buf,sizeof(T),nT,fp);
        free(buf);
    }
    fclose(fp);
}

template <class T>
bool BitSlicedIndex<T>::find(T *value, size_t *pos) {
    // The BitSlicedIndex holds sorted distinct values (for kmerizer)
    // initially value could be anywhere in the index
    BitVector<T> *range = new BitVector<T>();
    if (*pos > 0) {
        range->appendFill0(*pos);
        range->appendFill1(nValues - *pos);
    }
    else
        range->appendFill1(nValues);
    // refine the range
    for(int w=0;w<nwords;w++) {
        for(int b=0; b<nbits; b++) {
            // fprintf(stderr,"find() w %i b %i bit %c nwords %zi\n",w,b,((value[w] >> (nbits - b - 1)) & 1)?'1':'0', bvec[w*nbits + b]->getNWords());
            if ((value[w] >> (nbits - b - 1)) & 1) // value has a 1 here
                range->andYes(bvec[w*nbits + b]);
            else
                range->andNot(bvec[w*nbits + b]);
            if (range->getCount() == 0) return false;
        }
    }
    range->firstActiveWord(); // this fixed a seg fault that was probably a symptom of some other bug...
    return range->nextSetBit(pos);
}

#endif