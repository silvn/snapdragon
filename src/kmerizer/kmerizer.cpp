#include "kmerizer.h"
#include <bitset>
#include <algorithm> // stl sort
#include <cstdlib> // qsort()
#include <cstdio>  // sprintf()
#include <cstring> // memcpy()
#include <sys/stat.h> // mkdir()
#include <sys/time.h> // gettimeofday()

Kmerizer::Kmerizer(const size_t klength,
                   const size_t nthreads,
                   const char * outdir,
                   const char   mode)
{
    this->k            = klength-4;
    this->outdir       = new char[strlen(outdir)]; strcpy(this->outdir, outdir);
    this->mode         = mode;
    this->nthreads     = nthreads;

    this->kmask       = 0xFFFFFFFFFFFFFFFFULL;
    this->shiftlastby = 62;
    this->lshift      = 0;
    this->rshift      = 0;
    if (k % 32 > 0) {
        lshift      = 2 * (k % 32);
        rshift      = 64 - lshift;
        kmask       = (1ULL << lshift) - 1;
        shiftlastby = lshift - 2;
    }
    this->nwords     = ((k-1)>>5)+1;
    this->kmerSize   = this->nwords * sizeof(kword_t);

    tp.size_controller().resize(nthreads);
}

int Kmerizer::allocate(const size_t maximem) {
    fprintf(stderr, "Kmerizer::allocate(%zi)", maximem);
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);

    memset(batches,      0, sizeof(uint32_t) * NBINS);
    memset(kmerLutS,     0, sizeof(uint32_t) * NBINS);
    memset(kmerBufTally, 0, sizeof(uint32_t) * NBINS);

    uint32_t capacity = maximem / kmerSize / NBINS;
    int slots = 1*500/nwords; // 256*8*500 = one megabyte for the ProducerConsumerQueues?
    for (size_t i = 0; i < NBINS; i++) {
        // if      (nwords == 1) PCQ[i] = new folly::ProducerConsumerQueue<kmer1_t>(slots); 
        // else if (nwords == 2) PCQ[i] = new folly::ProducerConsumerQueue<kmer2_t>(slots);
        // else if (nwords == 3) PCQ[i] = new folly::ProducerConsumerQueue<kmer3_t>(slots);
        // else if (nwords == 4) PCQ[i] = new folly::ProducerConsumerQueue<kmer4_t>(slots);
        // else if (nwords == 5) PCQ[i] = new folly::ProducerConsumerQueue<kmer5_t>(slots);
        // else if (nwords == 6) PCQ[i] = new folly::ProducerConsumerQueue<kmer6_t>(slots);
        // else if (nwords == 7) PCQ[i] = new folly::ProducerConsumerQueue<kmer7_t>(slots);
        // else if (nwords == 8) PCQ[i] = new folly::ProducerConsumerQueue<kmer8_t>(slots);
        if      (nwords == 1) PCQ[i] = new boost::lockfree::spsc_queue<kmer1_t>(slots); 
        else if (nwords == 2) PCQ[i] = new boost::lockfree::spsc_queue<kmer2_t>(slots);
        else if (nwords == 3) PCQ[i] = new boost::lockfree::spsc_queue<kmer3_t>(slots);
        else if (nwords == 4) PCQ[i] = new boost::lockfree::spsc_queue<kmer4_t>(slots);
        else if (nwords == 5) PCQ[i] = new boost::lockfree::spsc_queue<kmer5_t>(slots);
        else if (nwords == 6) PCQ[i] = new boost::lockfree::spsc_queue<kmer6_t>(slots);
        else if (nwords == 7) PCQ[i] = new boost::lockfree::spsc_queue<kmer7_t>(slots);
        else if (nwords == 8) PCQ[i] = new boost::lockfree::spsc_queue<kmer8_t>(slots);

        kmerBufSize[i] = capacity;
        kmerBuf[i] = (kword_t *) calloc(capacity, kmerSize);
        if (kmerBuf[i] == NULL) return 1;

        // heuristic: expect the A* kmer to be common
        kmerLutS[i] = 1;
        kmerLutK[i] = (kword_t *) calloc(1, kmerSize);
        if (kmerLutK[i] == NULL) return 1;
        kmerLutV[i] = (uint32_t *) calloc(1, sizeof(uint32_t));
        if (kmerLutV[i] == NULL) return 1;
    }
    state = COUNT;
    // schedule the workers
    for(size_t i=0; i<nthreads; i++)
             if (nwords == 1) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer1_t>, this, i ) );
        else if (nwords == 2) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer2_t>, this, i ) );
        else if (nwords == 3) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer3_t>, this, i ) );
        else if (nwords == 4) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer4_t>, this, i ) );
        else if (nwords == 5) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer5_t>, this, i ) );
        else if (nwords == 6) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer6_t>, this, i ) );
        else if (nwords == 7) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer7_t>, this, i ) );
        else if (nwords == 8) tp.schedule( boost::bind( &Kmerizer::consumeKmers<kmer8_t>, this, i ) );

    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    fprintf(stderr," took %f seconds\n",elapsedTime);
    return 0;
}

bool operator<(const kmer1_t& a, const kmer1_t& b) {
    return a.a < b.a;
}
bool operator<(const kmer2_t& a, const kmer2_t& b) {
    int i=0, s=2;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
}
bool operator<(const kmer3_t& a, const kmer3_t& b) {
    int i=0, s=3;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
}
bool operator<(const kmer4_t& a, const kmer4_t& b) {
    int i=0, s=4;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
}
bool operator<(const kmer5_t& a, const kmer5_t& b) {
    int i=0, s=5;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
}
bool operator<(const kmer6_t& a, const kmer6_t& b) {
    int i=0, s=6;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
}
bool operator<(const kmer7_t& a, const kmer7_t& b) {
    int i=0, s=7;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
}
bool operator<(const kmer8_t& a, const kmer8_t& b) {
    int i=0, s=8;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
    // for (size_t i=0;i<8;i++) {
    //     if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
    //     if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    // }
    // return false;
}

template<class T>
kword_t * Kmerizer::canonicalize(kword_t *packed, kword_t *rcpack, size_t *bin) const {

    // copy the words (reverse their order)
    for (size_t i=0;i<nwords;i++)
        rcpack[i] = packed[nwords-1-i];

    // shift bits from word to word if necessary
    if (lshift) {
        for (size_t i=0;i<nwords-1;i++) {
            rcpack[i] |= rcpack[i+1] << lshift;
            rcpack[i+1] >>= rshift;
        }
    }

    size_t rcbin = revcomp_8(rcpack[0]);
    if (rcbin > *bin)
        return packed;

    // shift each word by 8 bits
    rcpack[0] >>= 8;
    for (size_t i=0;i<nwords-1;i++) {
        rcpack[i] |= rcpack[i+1] << 56;
        rcpack[i+1] >>= 8;
    }
    if (lshift >= 8)
        rcpack[nwords-1] |= *bin << (lshift-8);
    else if (lshift > 0) {
        rcpack[nwords-2] |= *bin << 56+lshift;
        rcpack[nwords-1] |= *bin >> lshift;
    }
    else
        rcpack[nwords-1] |= *bin << 56;

    int cmp=0;
    if (rcbin > *bin)
        cmp = 1;
    for (size_t i=0;i<nwords;i++) {
        rcpack[i] = revcomp(rcpack[i]);
        if (i==nwords-1) rcpack[i] >>=rshift;
        if (cmp == 0) {
            if (packed[i] < rcpack[i])
                cmp = -1;
            else if (packed[i] > rcpack[i])
                cmp = 1;
        }
    }
    if(cmp > 0) {
        *bin = rcbin;
        return rcpack;
    }
    else
        return packed;
}


// version for nwords==1

kword_t * Kmerizer::canonicalize1(kword_t *packed, kword_t *rcpack, size_t *bin) const {

    *rcpack = *packed;

    if (k < 4)
        // not enough bits at the end of the word to fill a bin
        // so prepend the bin's 8 bits now
        *rcpack |= *bin << lshift;
    
    size_t rcbin = revcomp_8(*rcpack);
    if (rcbin > *bin)
        return packed;

    // shift word by 8 bits
    *rcpack >>= 8;

    // prepend the bin's 8 bits
    if (k >= 4) {
        if (lshift == 0)
            *rcpack |= *bin << 56;
        else
            *rcpack |= *bin << (lshift - 8);
    }
    *rcpack = revcomp(*rcpack);
    *rcpack >>=rshift;

    if (rcbin > *bin) {
        *bin = rcbin;
        return rcpack;
    }
    if (*packed < *rcpack)
        return packed;
    else if (*packed > *rcpack)
        return rcpack;
}


// nwords > 1
template<class T>
size_t Kmerizer::nextKmer(kword_t* kmer, size_t bin, const char nucl) {
    // update the bin
    bin <<= 2;
    bin |= (kmer[0] >> 62);
    bin &= 255; // mask out the higher order bits
    
    // shift first kmer by 2 bits
    kmer[0] <<= 2;
    for (size_t w = 1; w < nwords - 1; w++) { // middle (full length) words
        kmer[w-1] |= (kmer[w] >> 62); // move the top two bits
        kmer[w] <<= 2; // make room for the next word
    }
    // last word (if not the first)
    kmer[nwords-2] |= (kmer[nwords-1] >> shiftlastby);
    kmer[nwords-1] <<= 2;
    kmer[nwords-1] |= twoBit(nucl);
    kmer[nwords-1] &= kmask;
    return bin;
}

// nwords == 1

size_t Kmerizer::nextKmer1(kword_t* kmer, size_t bin, const char nucl) {
    // update the bin
    bin <<= 2;
    bin |= (*kmer >> shiftlastby);
    bin &= 255; // mask out the higher order bits
    
    // shift first kmer by 2 bits
    *kmer <<= 2;
    *kmer |= twoBit(nucl);
    *kmer &= kmask;
    return bin;
}


void Kmerizer::addSequence(const char* seq, const int length) {
    if      (nwords == 1) addSeq1 (seq, length);
    else if (nwords == 2) addSeq<kmer2_t> (seq, length);
    else if (nwords == 3) addSeq<kmer3_t> (seq, length);
    else if (nwords == 4) addSeq<kmer4_t> (seq, length);
    else if (nwords == 5) addSeq<kmer5_t> (seq, length);
    else if (nwords == 6) addSeq<kmer6_t> (seq, length);
    else if (nwords == 7) addSeq<kmer7_t> (seq, length);
    else if (nwords == 8) addSeq<kmer8_t> (seq, length);
}

template<class T>
void Kmerizer::addSeq(const char* seq, const int length) {
    kword_t packed[nwords];
    kword_t rcpack[nwords];
    memset(packed, 0, kmerSize);
    size_t bin=0;
    for (size_t i = 0; i < k+4; i++)
        bin = nextKmer<T> (packed,bin,seq[i]);

    kword_t *kmer = packed;
    if (mode == CANONICAL)
        kmer = canonicalize<T> (packed,rcpack,&bin);

    // while (! ((folly::ProducerConsumerQueue<T>*)PCQ[bin])->write(*((T*)kmer))) continue;
    while (! ((boost::lockfree::spsc_queue<T>*)PCQ[bin])->push(*((T*)kmer))) ;
    // pack the rest of the sequence
    for (int i=k+4; i<length;i++) {
        bin = nextKmer<T> (packed, bin, seq[i]);
        if (mode == CANONICAL)
            kmer = canonicalize<T> (packed,rcpack, &bin);
        // insertKmer(bin, kmer);
        // while (! ((folly::ProducerConsumerQueue<T>*)PCQ[bin])->write(*((T*)kmer))) continue;
        while (! ((boost::lockfree::spsc_queue<T>*)PCQ[bin])->push(*((T*)kmer))) ;
    }
}

void Kmerizer::addSeq1(const char* seq, const int length) {
    kword_t packed=0;
    kword_t rcpack=0;
    size_t bin=0;
    for (size_t i = 0; i < k+4; i++)
        bin = nextKmer1(&packed,bin,seq[i]);

    kword_t *kmer = &packed;
    if (mode == CANONICAL)
        kmer = canonicalize1(&packed,&rcpack,&bin);

    // while (! ((folly::ProducerConsumerQueue<kmer1_t>*)PCQ[bin])->write(*((kmer1_t*)kmer))) continue;
    while (! ((boost::lockfree::spsc_queue<kmer1_t>*)PCQ[bin])->push(*((kmer1_t*)kmer))) ;

    // pack the rest of the sequence
    for (int i=k+4; i<length;i++) {
        bin = nextKmer1(&packed, bin, seq[i]);
        if (mode == CANONICAL)
            kmer = canonicalize1(&packed, &rcpack, &bin);
        // while (! ((folly::ProducerConsumerQueue<kmer1_t>*)PCQ[bin])->write(*((kmer1_t*)kmer))) continue;
        while (! ((boost::lockfree::spsc_queue<kmer1_t>*)PCQ[bin])->push(*((kmer1_t*)kmer))) ;
    }
}


template<class T>
void Kmerizer::consumeKmers(size_t threadId) {
    T kmer;
    bool empty;
    size_t bin;
    size_t attempts=0;
    size_t hits = 0;
    do {
        empty = true;
        bin = threadId;
        while (bin < NBINS) {
            attempts++;
            if (((boost::lockfree::spsc_queue<T>*)PCQ[bin])->pop(kmer)) {
                hits++;
                empty = false;
                if (nwords == 1)
                    insertKmer1(bin, (kword_t *) &kmer);
                else
                    insertKmer(bin, (kword_t *) &kmer);
            }
            else 
                bin+=nthreads;
        }
    } while (state == COUNT || !empty);
    fprintf(stderr,"consumeKmers(%zi) attempts: %zi hits: %zi\n",threadId, attempts, hits);
}
// template<class T>
// void Kmerizer::consumeKmers(size_t threadId) {
//     T* kmer;
//     bool empty;
//     size_t bin;
//     size_t attempts=0;
//     size_t hits = 0;
//     do {
//         empty = true;
//         bin = threadId;
//         while (bin < NBINS) {
//             kmer = ((folly::ProducerConsumerQueue<T>*)PCQ[bin])->frontPtr();
//             attempts++;
//             if (kmer) {
//                 hits++;
//                 empty = false;
//                 if (nwords == 1)
//                     insertKmer1(bin, (kword_t *) kmer);
//                 else
//                     insertKmer(bin, (kword_t *) kmer);
//                 ((folly::ProducerConsumerQueue<T>*)PCQ[bin])->popFront();
//             }
//             else 
//                 bin+=nthreads;
//         }
//     } while (state == COUNT || !empty);
//     fprintf(stderr,"consumeKmers(%zi) attempts: %zi hits: %zi\n",threadId, attempts, hits);
// }

void Kmerizer::insertKmer(size_t bin, kword_t *kmer) {

    int rc = searchLut(kmer,bin);
    if (rc >= 0) // found!
        kmerLutV[bin][rc]++;
    
    else { // kmer is not in the lookup table
        uint32_t offset = kmerBufTally[bin];
        if (offset < kmerBufSize[bin]) {
            // memcpy the kmer to this location
            memcpy(kmerBuf[bin] + offset*nwords, kmer, kmerSize); // is specialization worth it for kmer1_t?
            kmerBufTally[bin]++;
            if (offset + 1 == kmerBufSize[bin]) { // buffer is full
                uniqify(bin);
            }
        }
        else {
            fprintf(stderr,"the impossible happened\n");
            exit(3);
        }
    }
}

void Kmerizer::insertKmer1(size_t bin, kword_t *kmer) {

    int rc = searchLut1(kmer,bin);
    if (rc >= 0) // found!
        kmerLutV[bin][rc]++;
    
    else { // kmer is not in the lookup table
        uint32_t offset = kmerBufTally[bin];
        if (offset < kmerBufSize[bin]) {
            kmerBuf[bin][offset] = *kmer;
            kmerBufTally[bin]++;
            if (offset + 1 == kmerBufSize[bin]) { // buffer is full
                uniqify1(bin);
            }
        }
        else {
            fprintf(stderr,"the impossible happened\n");
            exit(3);
        }
    }
}

int kmercmp(const void *k1, const void *k2, size_t nwords) {
    for (size_t i=0;i<nwords;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}

int Kmerizer::searchLut(kword_t *kmer, size_t bin) {
    kword_t* A = kmerLutK[bin];
    int imin = 0;
    int imax = kmerLutS[bin];
    while (imin < imax) {
        int imid = (imin+imax)>>1;
        int cmp = kmercmp(A + nwords*imid, kmer, nwords);
//        int cmp = memcmp(A+nwords*imid, kmer, kmerSize);
        if (cmp < 0)
            imin = imid + 1;
        else if (cmp == 0)
            return imid;
        else
            imax = imid;
    }
    return -1;
}


int Kmerizer::searchLut1(kword_t *kmer, size_t bin) {
    kword_t* A = kmerLutK[bin];
    int imin = 0;
    int imax = kmerLutS[bin];
    while (imin < imax) {
        int imid = (imin+imax)>>1;
        if (A[imid] < *kmer)
            imin = imid + 1;
        else if(A[imid] == *kmer)
            return imid;
        else
            imax = imid;
    }
    return -1;
}



void Kmerizer::uniqify(size_t bin) {
//    fprintf(stderr,"uniqify(%zi) kmerLutS: %u, kmerBufTally: %u kmerBufSize: %u\n",bin,kmerLutS[bin],kmerBufTally[bin], kmerBufSize[bin]);
    if (kmerBufTally[bin] == 0) return;
    // sort the kmers in the buffer
    if      (nwords == 1) sort((kmer1_t*)kmerBuf[bin], (kmer1_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    else if (nwords == 2) sort((kmer2_t*)kmerBuf[bin], (kmer2_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    else if (nwords == 3) sort((kmer3_t*)kmerBuf[bin], (kmer3_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    else if (nwords == 4) sort((kmer4_t*)kmerBuf[bin], (kmer4_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    else if (nwords == 5) sort((kmer5_t*)kmerBuf[bin], (kmer5_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    else if (nwords == 6) sort((kmer6_t*)kmerBuf[bin], (kmer6_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    else if (nwords == 7) sort((kmer7_t*)kmerBuf[bin], (kmer7_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    else if (nwords == 8) sort((kmer8_t*)kmerBuf[bin], (kmer8_t*)(kmerBuf[bin] + kmerBufTally[bin]));
    
    // uniq -c
    uint32_t distinct = 0;
    vector<uint32_t> tally;
    tally.push_back(1); // first kmer
    for (size_t i=1;i<kmerBufTally[bin];i++) {
        kword_t *ith = kmerBuf[bin] + i*nwords;
        if(kmercmp(kmerBuf[bin] + distinct*nwords, ith, nwords) == 0)
        // if (memcmp(kmerBuf[bin] + distinct*nwords, ith, kmerSize) == 0)
            tally.back()++;
        else {
            distinct++;
            tally.push_back(1);
            memcpy(kmerBuf[bin] + distinct*nwords, ith, kmerSize);
        }
    }
    kmerBufTally[bin] = distinct+1;
    // update lookup table - this realloc's the kmer lut and buffer
    if (kmerLutS[bin] == 1)
        updateLut(bin,tally);
    // output distinct kmers
    if (kmerBufTally[bin] > 0) {
        batches[bin]++;
        // writeBatch(bin);
        // declare the buffer empty
        kmerBufTally[bin]=0;
    }
}

void Kmerizer::uniqify1(size_t bin) {
//    fprintf(stderr,"uniqify(%zi) kmerLutS: %u, kmerBufTally: %u kmerBufSize: %u\n",bin,kmerLutS[bin],kmerBufTally[bin], kmerBufSize[bin]);
    if (kmerBufTally[bin] == 0) return;
    // sort the kmers in the buffer
    sort(kmerBuf[bin], (kmerBuf[bin] + kmerBufTally[bin]));
    
    // uniq -c
    uint32_t distinct = 0;
    vector<uint32_t> tally;
    tally.push_back(1); // first kmer
    for (size_t i=1;i<kmerBufTally[bin];i++) {
        if(kmerBuf[bin][distinct] == kmerBuf[bin][i])
        // if (memcmp(kmerBuf[bin] + distinct*nwords, ith, kmerSize) == 0)
            tally.back()++;
        else {
            distinct++;
            tally.push_back(1);
            kmerBuf[bin][distinct] = kmerBuf[bin][i];
        }
    }
    kmerBufTally[bin] = distinct+1;
    // update lookup table - this realloc's the kmer lut and buffer
    if (kmerLutS[bin] == 1)
        updateLut(bin,tally);
    // output distinct kmers
    if (kmerBufTally[bin] > 0) {
        batches[bin]++;
        // writeBatch(bin);
        // declare the buffer empty
        kmerBufTally[bin]=0;
    }
}

void Kmerizer::writeBatch(size_t bin, vector<uint32_t> &tally) {
    // assume writing to disk is faster than compressing and writing to disk
    //
    // there is always a trade-off with "big data" when reading/writing to disk
    // if the IO is fast, then it doesn't make sense to spend CPU cycles
    // compressing and decompressing the data
    // merging batches of kmers can be done independantly, and doesn't use
    // as much memory as counting kmers since it is done by iterating over the kmers
    // However, if you are using a cluster to merge all 256 bins in parallel
    // the IO speed over the network could be a limiting factor.
    // Therefore, we should offer the choice to the user of whether to compress
    // each batch or write the raw data and compress while merging.
    FILE *fp;
    char kmer_file[100];
    sprintf(kmer_file,"%s/%zi-mers.%zi.%u.raw",outdir,k,bin,batches[bin]);
    fp = fopen(kmer_file, "wb");
    fwrite(kmerBuf[bin], kmerSize, kmerBufTally[bin], fp);
    fclose(fp);
    // convert tally to a range encoded bitmap self-index?
    // write to disk
    char kmerCount_file[100];
    sprintf(kmerCount_file,"%s/%zi-mers.%zi.%u.count",outdir,k,bin,batches[bin]);
    fp = fopen(kmerCount_file, "wb");
    fwrite(tally.data(), sizeof(uint32_t), kmerBufTally[bin], fp);
    fclose(fp); 
}
void Kmerizer::updateLut(size_t bin, vector<uint32_t> &tally) {
    // copy everything into the lookup table
    kmerLutK[bin] = (kword_t*) realloc(kmerLutK[bin], (kmerLutS[bin] + tally.size())*kmerSize);
    kmerLutV[bin] = (uint32_t*) realloc(kmerLutV[bin], (kmerLutS[bin] + tally.size())*sizeof(uint32_t));
    memcpy(kmerLutK[bin] + kmerLutS[bin]*nwords, kmerBuf[bin], tally.size()*kmerSize);
    memcpy(kmerLutV[bin] + kmerLutS[bin], tally.data(), tally.size()*sizeof(uint32_t));
    kmerLutS[bin] += tally.size();
    kmerBufSize[bin] -= tally.size();
    kmerBuf[bin] = (kword_t*) realloc(kmerBuf[bin], kmerBufSize[bin]*kmerSize);
    kmerBufTally[bin] = 0;
}


/*
void Kmerizer::updateLut(size_t bin, vector<uint32_t> &tally) {
    uint32_t total=0;
    uint32_t oversum=0;
    vector<uint32_t> counts (256,0);
    vector<uint32_t>::iterator it;
    for (it = tally.begin(); it < tally.end(); ++it) {
        total += *it;
        if (*it < 256)
            counts[*it]++;
        else
            oversum += *it;
    }
    total >>= 1; // half of the total
    uint32_t cutoff=256;
    while (oversum < total && cutoff > 1) {
        cutoff--;
        oversum += counts[cutoff]*cutoff;
    }
    if (cutoff>0) {
        // save the positions of common kmers in the buffer (and tally)
        vector<uint32_t> common;
        uint32_t common_count=0;
        for(uint32_t i=0;i<tally.size();i++)
            if (tally[i] > cutoff) {
                common_count++;
                if (common.size()==0 || i > common.back()+1) {
                    common.push_back(i);
                    common.push_back(i);
                }
                else
                    common.back()++;
            }
        assert(common_count>0);
        // reallocate space in the lookup table
        kmerLutK[bin] = (kword_t*) realloc(kmerLutK[bin], (common_count + kmerLutS[bin])*kmerSize);
        if (kmerLutK[bin] == NULL) {
            // error in realloc
            fprintf(stderr,"error in kmerLutK realloc");
            exit(3);
        }
        kmerLutV[bin] = (uint32_t*) realloc(kmerLutV[bin], (common_count + kmerLutS[bin])*sizeof(uint32_t));
        if (kmerLutV[bin] == NULL) {
            // error in realloc
            fprintf(stderr,"error in kmerLutV realloc");
            exit(3);
        }
        // memcpy frequent kmers and their counts into the lookup table
        // and shift rare kmers and their counts to fill the void
        uint32_t *tarray = tally.data();
        // copy the first run of common kmers into the lookup table
        uint32_t runlen = common[1] - common[0] + 1;
        memcpy(kmerLutK[bin] + kmerLutS[bin]*nwords, kmerBuf[bin] + common[0]*nwords, runlen*kmerSize);
        memcpy(kmerLutV[bin] + kmerLutS[bin], tarray + common[0], runlen*sizeof(uint32_t));
        kmerLutS[bin] += runlen;
        uint32_t n_rare = common[0];
        for(uint32_t i=2;i<common.size();i+=2) {
            // fill the void left by copying the previous run of common kmers to the lookup table
            runlen = common[i] - common[i-1] - 1;
            memmove(kmerBuf[bin] + n_rare*nwords, kmerBuf[bin] + (common[i-1] + 1)*nwords, runlen*kmerSize);
            memmove(tarray + n_rare, tarray + common[i-1]+1, runlen*sizeof(uint32_t));
            n_rare += runlen;
            // copy this run of common kmers to the lookup table
            runlen = common[i+1] - common[i] + 1;
            memcpy(kmerLutK[bin] + kmerLutS[bin]*nwords, kmerBuf[bin] + common[i]*nwords, runlen*kmerSize);
            memcpy(kmerLutV[bin] + kmerLutS[bin], tarray + common[i], runlen*sizeof(uint32_t));
            kmerLutS[bin] += runlen;
        }
        // fill the void left by the last run of common kmers to the lookup table
        runlen = tally.size() - common.back() - 1;
        if (runlen > 0) { // in case the last kmer was moved
            memmove(kmerBuf[bin] + n_rare*nwords, kmerBuf[bin] + (common.back() + 1)*nwords, runlen*kmerSize);
            memmove(tarray + n_rare, tarray + common.back()+1, runlen*sizeof(uint32_t));
            n_rare += runlen;
        }
        // reallocate the kmer buffer and the corresponding tally vector
        kmerBuf[bin] = (kword_t*) realloc(kmerBuf[bin], n_rare*kmerSize);
        kmerBufSize[bin] = n_rare;
        kmerBufTally[bin] = n_rare;
        tally.resize(n_rare);
    }
}
*/

void Kmerizer::save() {
    state = SAVE;
    tp.wait();
    fprintf(stderr,"counting finished, scheduling saveBin()\n");
    for(size_t bin=0; bin < NBINS; bin++)
        tp.schedule( boost::bind( &Kmerizer::saveBin, this, bin ) );

    // flush the task queue
    tp.wait();
    fprintf(stderr,"all bins saved to disk\n");
    state = MERGE;
    // merge the batches
    for(size_t bin=0; bin < NBINS; bin++)
        tp.schedule( boost::bind( &Kmerizer::mergeBin, this, bin ) );
    tp.wait();
    fprintf(stderr,"all bins merged\n");
    state = QUERY;
}

void Kmerizer::saveBin(size_t bin) {
    // all we want to do here is run uniqify and write the lookup table to disk.
    // merging is done separately because bins can be merged independantly
    if (nwords == 1) uniqify1(bin);
    else  uniqify(bin);
//    return;
    // write the common kmers lut to disk
    if (0) { // kmerLutS[bin] > 1 || kmerLutV[bin][0]>0) 
        FILE *fp;
        char kmer_file[100];
        sprintf(kmer_file,"%s/%zi-mers.%zi.0.raw",outdir,k,bin);
        fp = fopen(kmer_file, "wb");
        fwrite(kmerLutK[bin], kmerSize, kmerLutS[bin], fp);
        fclose(fp);
        // convert tally to a range encoded bitmap self-index?
        char kmerCount_file[100];
        sprintf(kmerCount_file,"%s/%zi-mers.%zi.0.count",outdir,k,bin);
        fp = fopen(kmerCount_file, "wb");
        fwrite(kmerLutV[bin], sizeof(uint32_t), kmerLutS[bin], fp);
        fclose(fp);        
        // declare the lookup table empty
        kmerLutS[bin]=0;
    }
}

void Kmerizer::mergeBin(size_t bin) {
    // create new bitmap indexes for the kmers (binary encoded) and their frequencies (range encoded)
    // open each batch of files (mmap)
    // review notes on hybrid merge
}