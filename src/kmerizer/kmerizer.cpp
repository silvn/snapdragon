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
    this->maxSeqLength = 100000;

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
    memset(usingLut, 0, sizeof(int) * NBINS);
    memset(usingBuf, 0, sizeof(int) * NBINS);

    uint32_t capacity = maximem / kmerSize / NBINS;
    for (size_t i = 0; i < NBINS; i++) {
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
    state = COUNTING;

    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    fprintf(stderr," took %f seconds\n",elapsedTime);
    return 0;
}

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
    int klength = k+4;
    if (nwords == 1)
        for(int offset=0; offset + klength <= length; offset += maxSeqLength - klength) {
            int sslength = (offset + maxSeqLength > length) ? length - offset : maxSeqLength;
//            pool->run_task( boost::bind( &Kmerizer::addSeq1, this, seq + offset, sslength ) );
            tp.schedule( boost::bind( &Kmerizer::addSeq1, this, seq + offset, sslength ) );
        }
    else
        for(int offset=0; offset + klength <= length; offset += maxSeqLength - klength) {
            int sslength = (offset + maxSeqLength > length) ? length - offset : maxSeqLength;
//            pool->run_task( boost::bind( &Kmerizer::addSeq, this, seq + offset, sslength ) );
            tp.schedule( boost::bind( &Kmerizer::addSeq, this, seq + offset, sslength ) );
        }
    tp.wait();
}

void Kmerizer::addSeq(const char* seq, const int length) {
    kword_t packed[nwords];
    kword_t rcpack[nwords];
    memset(packed, 0, kmerSize);
    size_t bin=0;
    for (size_t i = 0; i < k+4; i++)
        bin = nextKmer(packed,bin,seq[i]);

    kword_t *kmer = packed;
    if (mode == CANONICAL)
        kmer = canonicalize(packed,rcpack,&bin);

    insertKmer(bin, kmer);

    // pack the rest of the sequence
    for (int i=k+4; i<length;i++) {
        bin = nextKmer(packed, bin, seq[i]);
        if (mode == CANONICAL)
            kmer = canonicalize(packed,rcpack, &bin);
        insertKmer(bin, kmer);
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

    insertKmer1(bin, kmer);

    // pack the rest of the sequence
    for (int i=k+4; i<length; i++) {
        bin = nextKmer1(&packed, bin, seq[i]);
        if (mode == CANONICAL)
            kmer = canonicalize1(&packed,&rcpack, &bin);
        insertKmer1(bin, kmer);
    }
}

void requestResource(int* semaphore) {
    int ncount = *semaphore, count;
    do {
        while (ncount < 0) {
            usleep(1);
            ncount = *semaphore;
        }
        count = ncount;
        ncount = __sync_val_compare_and_swap(semaphore, count, count + 1);
    } while (ncount != count);
}

void releaseResource(int* semaphore) {
    int rc = __sync_sub_and_fetch(semaphore,1);
    assert(rc >= 0);
}

void Kmerizer::insertKmer(size_t bin, kword_t *kmer) {
    requestResource(usingLut + bin);
    int rc = searchLut(kmer,bin);
    if (rc >= 0) { // found!
        uint32_t x = __sync_add_and_fetch(kmerLutV[bin]+rc,1);
        releaseResource(usingLut + bin);
    }
    else { // kmer is not in the lookup table
        releaseResource(usingLut + bin);
        uint32_t offset = __sync_fetch_and_add(kmerBufTally + bin,1);
        if (offset < kmerBufSize[bin]) {
            requestResource(usingBuf + bin);
            // memcpy the kmer to this location
            memcpy(kmerBuf[bin] + offset*nwords, kmer, kmerSize);
            releaseResource(usingBuf + bin);
            if (offset + 1 == kmerBufSize[bin]) { // buffer is full
                // lock the lookup table and the buffer by setting the semaphores to -1
                // but wait until the other threads have finished
                while (!__sync_bool_compare_and_swap(usingLut + bin, 0, -1)) { usleep(1); };
                while (!__sync_bool_compare_and_swap(usingBuf + bin, 0, -1)) { usleep(1); };
                // uniqify the bin
                uniqify(bin);
                // unlock the kmer buffer and lookup table (in that order)
                usingBuf[bin]=0;
                usingLut[bin]=0;
            }
        }
        else {
            // another thread has already taken the last slot
            // give up an try again later
            usleep(1); // 1 millisecond?
            insertKmer(bin, kmer);
        }
    }
}
void Kmerizer::insertKmer1(size_t bin, kword_t *kmer) {
    insertKmer(bin,kmer);
    // when insertKmer() is working properly, make a nwords==1 version
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


struct kmer1_t { kword_t first_word; };
struct kmer2_t { kword_t first_word; };
struct kmer3_t { kword_t first_word; };
struct kmer4_t { kword_t first_word; };
struct kmer5_t { kword_t first_word; };
struct kmer6_t { kword_t first_word; };
struct kmer7_t { kword_t first_word; };
struct kmer8_t { kword_t first_word; };
bool operator<(const kmer1_t& a, const kmer1_t& b) {
    return a.first_word < b.first_word;
}
bool operator<(const kmer2_t& a, const kmer2_t& b) {
    for (size_t i=0;i<2;i++) {
        if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
        if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    }
    return false;
}
bool operator<(const kmer3_t& a, const kmer3_t& b) {
    for (size_t i=0;i<3;i++) {
        if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
        if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    }
    return false;
}
bool operator<(const kmer4_t& a, const kmer4_t& b) {
    for (size_t i=0;i<4;i++) {
        if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
        if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    }
    return false;
}
bool operator<(const kmer5_t& a, const kmer5_t& b) {
    for (size_t i=0;i<5;i++) {
        if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
        if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    }
    return false;
}
bool operator<(const kmer6_t& a, const kmer6_t& b) {
    for (size_t i=0;i<6;i++) {
        if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
        if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    }
    return false;
}
bool operator<(const kmer7_t& a, const kmer7_t& b) {
    for (size_t i=0;i<7;i++) {
        if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
        if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    }
    return false;
}
bool operator<(const kmer8_t& a, const kmer8_t& b) {
    for (size_t i=0;i<8;i++) {
        if (*(&a.first_word+i) < *(&b.first_word+i)) return true;
        if (*(&a.first_word+i) > *(&b.first_word+i)) return false;
    }
    return false;
}

void Kmerizer::uniqify(size_t bin) {
//    fprintf(stderr,"uniqify(%zi) kmerLutS: %u, kmerBufTally: %u kmerBufSize: %u\n",bin,kmerLutS[bin],kmerBufTally[bin], kmerBufSize[bin]);
    // sort the kmers in the buffer
    if (nwords == 1)      sort((kmer1_t*)kmerBuf[bin], (kmer1_t*)(kmerBuf[bin] + kmerBufTally[bin]));
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
            tally.back()++;
        else {
            distinct++;
            tally.push_back(1);
            memcpy(kmerBuf[bin] + distinct*nwords, ith, kmerSize);
        }
    }
    kmerBufTally[bin] = distinct+1;
//    if (state == SAVING) return;
    // update lookup table - this realloc's the kmer lut and buffer
    if (kmerLutS[bin] == 1)
        updateLut(bin,tally);
    // output distinct kmers
    if (kmerBufTally[bin] > 0) {
        batches[bin]++;

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
        if (1) {
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
        // declare the buffer empty
        kmerBufTally[bin]=0;
    }
}

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

void Kmerizer::save() {
    tp.wait();
    fprintf(stderr,"counting finished, scheduling saveBin()\n");
//    state = SAVING;
    for(size_t bin=0; bin < NBINS; bin++)
//        pool->run_task( boost::bind( &Kmerizer::saveBin, this, bin ) );
        tp.schedule( boost::bind( &Kmerizer::saveBin, this, bin ) );

    // flush the task queue
    tp.wait();
    fprintf(stderr,"all saves finished\n");
}

void Kmerizer::saveBin(size_t bin) {
    // all we want to do here is run uniqify and write the lookup table to disk.
    // merging is done separately because bins can be merged independantly
    uniqify(bin);
//    return;
    if (kmerLutS[bin] > 1 || kmerLutV[bin][0]>0) {
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
