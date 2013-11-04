// #include "bitmap.h"


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
    // fread() the number of words
    FILE *fp;
    fp = fopen(fname, "rb");
    size_t result = fread(&nwords,sizeof(int),1,fp);
    if (result != 1) {fputs ("Reading error",stderr); exit (3);}
    result = fread(&nValues, sizeof(size_t),1,fp);
    if (result != 1) {fputs ("Reading error",stderr); exit (3);}
    int nVectors = nwords*nbits;
    bvec = (BitVector<T>**) malloc(nVectors*sizeof(BitVector<T>*));
    for(int i=0;i<nVectors;i++) {
        T nT,*bvbuffer;
        result = fread(&nT,sizeof(size_t),1,fp);
        if (result != 1) {fputs ("Reading error",stderr); exit (3);}
        bvbuffer = (T*) malloc(nT*sizeof(T));
        result = fread(bvbuffer,sizeof(T), nT, fp);
        bvec[i] = new BitVector<T>(bvbuffer);
        free(bvbuffer);
    }
    fclose(fp);
    // fprintf(stderr,"BitSlicedIndex(%s) nwords: %i, nValues: %zi\n",fname,nwords,nValues);

    this->bufferCapacity = nbits;
    this->bufferOffset=0;
    this->bufferStart=0;

    this->buffer = (T**) malloc(nwords*sizeof(T*));
    for(int i=0;i<nwords;i++)
        buffer[i] = (T*) malloc(bufferCapacity*sizeof(T));
    fillBuffer(0);
}

template <class T>
void BitSlicedIndex<T>::transpose(uint64_t A[64]) {
    int j, k;
    uint64_t m, t;
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
void BitSlicedIndex<T>::transpose(uint16_t A[16]) {
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
void BitSlicedIndex<T>::append(T* value) {
    for(int i=0; i<nwords; i++)
        buffer[i][bufferOffset] = value[i];
    bufferOffset++;
    if (bufferOffset == nbits) {
        for(int i=0;i<nwords;i++) {
            this->transpose(buffer[i]);
            for(int j=0;j<nbits;j++)
                bvec[i*nbits + j]->appendWord(buffer[i][j]);
        }
        bufferOffset = 0;
    }
    nValues++;
}

template <class T>
void BitSlicedIndex<T>::append(T* value, int x) {
    for(int i=0; i<nwords; i++)
        buffer[i][bufferOffset] = value[i];
    bufferOffset++;
    if (bufferOffset == nbits) {
        for(int i=0;i<nwords;i++) {
            // for(int j=0;j<nbits;j++) {
            //     fprintf(stderr,"%zi %i %i ",bufferOffset,i,j);
            //     printBits(sizeof(T),buffer[i]+j);
            //     fprintf(stderr,"\n");
            // }
            this->transpose(buffer[i]);
            // for(int j=0;j<nbits;j++) {
            //     fprintf(stderr,"%zi %i %i ",bufferOffset,i,j);
            //     printBits(sizeof(T),buffer[i]+j);
            //     fprintf(stderr,"\n");
            // }
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
    bufferStart = idx & ~(nbits - 1);
    for(int i=0;i<nwords;i++) {
        for(int j=0;j<nbits;j++) {
            bvec[i*nbits+j]->inflateWord(buffer[i]+j,bufferStart);
        }
        
        // for(int j=0;j<nbits;j++) {
        //     printBits(sizeof(T),buffer[i]+j);
        //     fprintf(stderr," %i\n",j);
        // }
        
        this->transpose(buffer[i]);

        // for(int j=0;j<nbits;j++) {
        //     printBits(sizeof(T),buffer[i]+j);
        //     fprintf(stderr," %i\n",j);
        // }
    }
}

template <class T>
bool BitSlicedIndex<T>::decode(size_t idx, T *value) {
    if (idx >= nValues) return false;
    if (idx >= (bufferStart + nbits) || idx < bufferStart)
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
            this->transpose(buffer[i]);
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
        // fprintf(stderr,"fwrite(buf,%i,%i,fp)\n",sizeof(T),nT);
        fwrite(buf,sizeof(T),nT,fp);
        // BitVector<T>* bv_test = new BitVector<T>(buf);
        ////////////////////////////////////////////////////////////////////////////////////compare bv_test with bvec[i]
        // fprintf(stderr,"bvec[%i] nValues %zi\n",i,nValues);
        // T A,B;
        // for(int o=0;o<nValues;o+=nbits) {
        //     // fprintf(stderr,"calling bv_test->inflateWord(&A,%i)\n",o);
        //     bv_test->inflateWord(&A,o);
        //     // fprintf(stderr,"calling bvec[%i]->inflateWord(&B,%i)\n",i,o);
        //     bvec[i]->inflateWord(&B,o);
        //     if (A != B) {
        //         fprintf(stderr,"ERROR %llu != %llu\n",A,B);
        //         exit(4);
        //     }
        //     // else
        //         // fprintf(stderr,"word[%i] %llu == %llu\n",o,A,B);
        // }
        free(buf);
    }
    fclose(fp);
}

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
                ranges.push_back(i*64 + ffs(bits) - 1); // ffs() defined in bitvector.h. overloads __builtin_ffs(bits)
                bits &= bits-1;
            }
        }
        ranges.insert(ranges.end(), overflow.begin(), overflow.end());
    }
    fprintf(stderr,"ranges.size() = %zi nValues: %zi\n",ranges.size(),nValues);
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
            lb = lower_bound(lb, ranges.end(), *it);
            idx = lb - ranges.begin();
            for(int i=prevIdx+1; i<=idx; i++) {
                // fprintf(stderr,"bvec[%i]->appendFill0(%u)\n",i,currentPos-runStart[i]);
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

    // bufferStart=0;
    // bufferSize = 1 << 6; // tweak me
    // buffer = (uint32_t*) malloc(bufferSize * sizeof(uint32_t));
    // fillBuffer(0);
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
    fprintf(stderr,"RangeEncodedIndex(%s) nVectors: %zi nValues: %zi\n",fname,nVectors,nValues);
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
    fwrite(&nVectors, sizeof(size_t),1,fp);
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
    }
    fclose(fp);
}

void RangeEncodedIndex::fillBuffer(size_t idx) {
    bufferStart = idx & ~(bufferSize - 1);
    // everything is >= ranges[0], so initialize the buffer to that value
    memcpy(buffer, buffer0, bufferSize*sizeof(uint32_t));
    // fprintf(stderr,"ranges[0]:%u",ranges[0]);
    // for(int j=0;j<bufferSize;j++) fprintf(stderr," %u",buffer[j]);
    for(int i=1;i<ranges.size();i++) {
        bvec[i]->fillSetBits(bufferStart,buffer,ranges[i]);
        // fprintf(stderr,"\nranges[%i]:%u",i,ranges[i]);
        // for(int j=0;j<bufferSize;j++) fprintf(stderr," %u",buffer[j]);
    }
    // exit(1);
}

bool RangeEncodedIndex::decode(size_t idx, uint32_t *value) {
    if (idx >= nValues) return false;
    if ((idx >= bufferStart + bufferSize) || idx < bufferStart)
        fillBuffer(idx);
    
    int offset= idx - bufferStart;
    *value = buffer[offset];
    // if (idx == 0) {
    //     fprintf(stderr,"rei->decode(%zi)\n",idx);
    //     for(size_t i=0;i<bufferSize;i++) fprintf(stderr," %u",buffer[i]);
    //     fprintf(stderr," buffer[%i] %u\n", offset, buffer[offset]);
    //     fprintf(stderr,"ranges.size(): %zi",ranges.size());
    //     for(size_t i=0;i<ranges.size();i++) fprintf(stderr, " %u",ranges[i]);
    //     fprintf(stderr,"\n");
    // }
    return true;
}