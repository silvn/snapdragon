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

    memset(batches,        0, sizeof(uint32_t) * NBINS);
    // memset(kmerLutS,       0, sizeof(uint32_t) * NBINS);
    if (nwords>1)
        memset(kmerBufTally,   0, sizeof(uint32_t) * NBINS);
    // memset(kmerBufChecked, 0, sizeof(uint32_t) * NBINS);

    uint32_t capacity = maximem / kmerSize / NBINS;
    for (size_t i = 0; i < NBINS; i++) {
        kmerBufSize[i] = capacity;
        if (nwords == 1)
            kmerBuf1[i].reserve(capacity);
        else
            kmerBuf[i] = (kword_t *) calloc(capacity, kmerSize);
        // if (kmerBuf[i] == NULL) return 1;
    }
    state = COUNT;

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

kword_t * Kmerizer::canonicalize(kword_t *packed, kword_t *rcpack, size_t *bin) {

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
        rcpack[i] = reverse_complement(rcpack[i]);
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
    return packed;
}


// version for nwords==1
kword_t * Kmerizer::canonicalize1(kword_t *packed, kword_t *rcpack, size_t *bin) {

    *rcpack = *packed;

    if (k < 4)
        // not enough bits at the end of the word to fill a bin
        // so prepend the bin's 8 bits now
        *rcpack |= *bin << lshift;
    
    // size_t rcbin = reverse_complement(*rcpack) >> 56;
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
    *rcpack = reverse_complement(*rcpack);
    *rcpack >>= rshift;

    if (rcbin > *bin) {
        *bin = rcbin;
        return rcpack;
    }
    if (*rcpack < *packed)
        return rcpack;
    return packed;
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
    
    // kmer[nwords-1] |= twoBit(nucl);
    // kmer[nwords-1] &= kmask;
    static const kword_t lut[4] = {0UL,1UL,3UL,2UL};
    kmer[nwords-1] = (kmer[nwords-1] | lut[(nucl >> 1) & 3]) & kmask;
    return bin;
}

// nwords == 1
size_t Kmerizer::nextKmer1(kword_t* kmer, size_t bin, const char nucl) {
    // update the bin
    bin = ((bin << 2) | (*kmer >> shiftlastby)) & 255;
    static const kword_t lut[4] = {0UL,1UL,3UL,2UL};
    *kmer = ((*kmer << 2) | lut[(nucl >> 1) & 3]) & kmask;
    return bin;
}


void Kmerizer::addSequence(const char* seq, const int length) {
    if (length < k+4) return; // skip sequences that are too short
    if (nwords == 1)
        addSeq1(seq, length);
    else
        addSeq(seq, length);
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
        bin = nextKmer (packed, bin, seq[i]);
        if (mode == CANONICAL)
            kmer = canonicalize(packed,rcpack, &bin);
        insertKmer(bin, kmer);
    }
}

void Kmerizer::addSeq1(const char* seq, const int length) {
    // do the first kmer
    kword_t packed=0;
    size_t bin=0;
    static const kword_t lut[4] = {0UL,1UL,3UL,2UL};
    for (size_t i = 0; i < k+4; i++) {
        // bin = nextKmer1(&packed,bin,seq[i]);
        bin = ((bin << 2) | (packed >> shiftlastby)) & 255;
        packed = ((packed << 2) | lut[(seq[i] >> 1) & 3]) & kmask;
    }
    // pack the rest of the sequence
    if (mode == CANONICAL) {
        kword_t rcpack=0;
        kword_t *kmer = canonicalize1(&packed,&rcpack,&bin);
        // insertKmer1(bin, kmer);
        kmerBuf1[bin].push_back(*kmer);
        if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
            for(size_t i=0;i<NBINS;i++)
                tp.schedule( boost::bind( &Kmerizer::uniqify1, this, i ) );
            tp.wait();
        }

        for (int i=k+4; i<length;i++) {
            // bin = nextKmer1(&packed, bin, seq[i]);
            bin = ((bin << 2) | (packed >> shiftlastby)) & 255;
            packed = ((packed << 2) | lut[(seq[i] >> 1) & 3]) & kmask;
            kmer = canonicalize1(&packed,&rcpack,&bin);
            // insertKmer1(bin, kmer);
            kmerBuf1[bin].push_back(*kmer);
            if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
                for(size_t i=0;i<NBINS;i++)
                    tp.schedule( boost::bind( &Kmerizer::uniqify1, this, i ) );
                tp.wait();
            }
        }
    }
    else {
        // insertKmer1(bin, &packed);
        kmerBuf1[bin].push_back(packed);
        if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
            for(size_t i=0;i<NBINS;i++)
                tp.schedule( boost::bind( &Kmerizer::uniqify1, this, i ) );
            tp.wait();
        }
        for (int i=k+4; i<length;i++) {
            // bin = nextKmer1(&packed, bin, seq[i]);
            bin = ((bin << 2) | (packed >> shiftlastby)) & 255;
            packed = ((packed << 2) | lut[(seq[i] >> 1) & 3]) & kmask;
            // insertKmer1(bin, &packed);
            kmerBuf1[bin].push_back(packed);
            if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
                for(size_t i=0;i<NBINS;i++)
                    tp.schedule( boost::bind( &Kmerizer::uniqify1, this, i ) );
                tp.wait();
            }
        }
    }
}


void Kmerizer::insertKmer(size_t bin, kword_t *kmer) {
    uint32_t offset = kmerBufTally[bin];
    // memcpy the kmer to this location
    memcpy(kmerBuf[bin] + offset*nwords, kmer, kmerSize);
    kmerBufTally[bin]++;
    if (kmerBufTally[bin] == kmerBufSize[bin]) { // buffer is full
        for(size_t i=0;i<NBINS;i++)
            tp.schedule( boost::bind( &Kmerizer::uniqify, this, i ) );
        tp.wait();
    }
}

void Kmerizer::insertKmer1(size_t bin, kword_t *kmer) {
    // uint32_t offset = kmerBufTally[bin];
    // kmerBuf[bin][offset] = *kmer;
    // kmerBufTally[bin]++;
    // if (kmerBufTally[bin] == kmerBufSize[bin]) {
    kmerBuf1[bin].push_back(*kmer);
    if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
        for(size_t i=0;i<NBINS;i++)
            tp.schedule( boost::bind( &Kmerizer::uniqify1, this, i ) );
        tp.wait();
    }
}

int kmercmp(const void *k1, const void *k2, size_t nwords) {
    for (size_t i=0;i<nwords;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}


void Kmerizer::uniqify(size_t bin) {
    if (state == SAVE || kmerBufTally[bin] > 0.9*(float)kmerBufSize[bin]) {
        vector<uint32_t> tally;
        kmerBufTally[bin] = sortCount(kmerBuf[bin],kmerBufTally[bin],tally);
        // output distinct kmers
        if (kmerBufTally[bin] > 0) {
            batches[bin]++;
            writeBatch(bin, tally);
            // declare the buffer empty
            kmerBufTally[bin]=0;
        }
    }
}

void Kmerizer::uniqify1(size_t bin) {
    // if (state == SAVE || kmerBufTally[bin] > 0.9*(float)kmerBufSize[bin]) {
    if (state == SAVE || kmerBuf1[bin].size() > 0.9*(float)kmerBufSize[bin]) {
        vector<uint32_t> tally;
        // kmerBufTally[bin] = binSortCount1(kmerBuf[bin],kmerBufTally[bin],tally);
        binSortCount1(kmerBuf1[bin],tally);

        // output distinct kmers
        // if (kmerBufTally[bin] > 0) {
        if (kmerBuf1[bin].size() > 0) {
            batches[bin]++;
            writeBatch(bin, tally);
            // declare the buffer empty
            // kmerBufTally[bin]=0;
            kmerBuf1[bin].clear();
        }
//        kmerBufChecked[bin]=0;
    }
}

uint32_t Kmerizer::sortCount(kword_t *kb, uint32_t kbt, vector<uint32_t> &tally) {
    if      (nwords == 1) sort((kmer1_t*)kb, (kmer1_t*)(kb + kbt));
    else if (nwords == 2) sort((kmer2_t*)kb, (kmer2_t*)(kb + kbt));
    else if (nwords == 3) sort((kmer3_t*)kb, (kmer3_t*)(kb + kbt));
    else if (nwords == 4) sort((kmer4_t*)kb, (kmer4_t*)(kb + kbt));
    else if (nwords == 5) sort((kmer5_t*)kb, (kmer5_t*)(kb + kbt));
    else if (nwords == 6) sort((kmer6_t*)kb, (kmer6_t*)(kb + kbt));
    else if (nwords == 7) sort((kmer7_t*)kb, (kmer7_t*)(kb + kbt));
    else if (nwords == 8) sort((kmer8_t*)kb, (kmer8_t*)(kb + kbt));
    // uniq -c
    uint32_t distinct = 0;
    tally.push_back(1); // first kmer
    for (size_t i=1;i<kbt;i++) {
        kword_t *ith = kb + i*nwords;
        if(kmercmp(kb + distinct*nwords, ith, nwords) == 0)
        // if (memcmp(kmerBuf[bin] + distinct*nwords, ith, kmerSize) == 0)
            tally.back()++;
        else {
            distinct++;
            tally.push_back(1);
            memcpy(kb + distinct*nwords, ith, kmerSize);
        }
    }
    return distinct+1;
}

void Kmerizer::sortCount1(vector<kword_t> &kb, vector<uint32_t> &tally) {
    sort(kb.begin(), kb.end());
    // uniq -c
    uint32_t distinct = 0;
    register kword_t prevKmer = kb[distinct]; // to be cache friendly
    tally.push_back(1); // first kmer
    for(vector<kword_t>::iterator it = kb.begin()+1; it < kb.end(); it++) {
        if(*it == prevKmer)
            tally.back()++;
        else {
            distinct++;
            tally.push_back(1);
            kb[distinct] = *it;
            prevKmer = *it;
        }
    }
    kb.resize(distinct+1);
}

uint32_t Kmerizer::sortCount1(kword_t *kb, uint32_t kbt, vector<uint32_t> &tally) {
    sort(kb, kb + kbt);
    // uniq -c
    uint32_t distinct = 0;
    kword_t prevKmer = kb[distinct]; // to be cache friendly
    tally.push_back(1); // first kmer
    for (size_t i=1;i<kbt;i++) {
        if(kb[i] == prevKmer)
            tally.back()++;
        else {
            distinct++;
            tally.push_back(1);
            kb[distinct] = kb[i];
            prevKmer = kb[i];
        }
    }
    return distinct+1;
}

void Kmerizer::binSortCount1(vector<kword_t> &kb, vector<uint32_t> &tally) {
    int hbits = 12;
    int hbins = 1 << hbits; // 4096
    uint32_t *histogram = (uint32_t *) calloc(hbins, sizeof(uint32_t));
    uint32_t distinct = 0;

    int shiftby = 2*k - hbits;
    if (shiftby <= 0) {
        // direct hashing
        for(size_t i=0;i<kb.size();i++)
            histogram[kb[i]]++;

        for(uint32_t b = 0; b < hbins; b++)
            if (histogram[b] > 0) {
                kb[distinct] = (kword_t) b;
                distinct++;
                tally.push_back(histogram[b]);
            }
    }
    else {
        vector<kword_t> keys [hbins];
        
        // count the number of kmers in each bin
        // for(size_t i=0;i<kb.size();i++)
        for(vector<kword_t>::iterator it = kb.begin()+1; it < kb.end(); it++)
            histogram[*it >> shiftby]++;
        
        // allocate space for each bin
        for(int i=0; i<hbins; i++)
            if (histogram[i] > 0)
                keys[i].reserve(histogram[i]);
        free(histogram);

        // loop over the kmers again and copy them into their bins
        // for(uint32_t i=0;i<kb.size();i++)
        for(vector<kword_t>::iterator it = kb.begin()+1; it < kb.end(); it++)
            keys[*it >> shiftby].push_back(*it);
        
        kb.clear();
        for(uint32_t b = 0; b < hbins; b++) {
            if (keys[b].size() > 0) {
                sortCount1(keys[b],tally);
                kb.insert(kb.end(),keys[b].begin(),keys[b].end());
            }
        }
    }
}

uint32_t Kmerizer::binSortCount1(kword_t *kb, uint32_t kbt, vector<uint32_t> &tally) {
    int hbits = 12;
    int hbins = 1 << hbits; // 4096
    uint32_t *histogram = (uint32_t *) calloc(hbins, sizeof(uint32_t));
    uint32_t distinct = 0;

    int shiftby = 2*k - hbits;
    if (shiftby <= 0) {
        // direct hashing
        for(size_t i=0;i<kbt;i++)
            histogram[kb[i]]++;

        for(uint32_t b = 0; b < hbins; b++)
            if (histogram[b] > 0) {
                kb[distinct] = (kword_t) b;
                distinct++;
                tally.push_back(histogram[b]);
            }
    }
    else {
        vector<kword_t> keys [hbins];
        
        // count the number of kmers in each bin
        for(size_t i=0;i<kbt;i++)
            histogram[kb[i] >> shiftby]++;
        
        // allocate space for each bin
        for(int i=0; i<hbins; i++)
            if (histogram[i] > 0)
                keys[i].reserve(histogram[i]);
        free(histogram);

        // loop over the kmers again and copy them into their bins
        for(uint32_t i=0;i<kbt;i++)
            keys[kb[i] >> shiftby].push_back(kb[i]);
        
        for(uint32_t b = 0; b < hbins; b++) {
            if (keys[b].size() > 0) {
                kword_t *K = keys[b].data();
                uint32_t bin_distinct = sortCount1(K,keys[b].size(),tally);
                memcpy(kb + distinct, K, sizeof(kword_t) * bin_distinct);
                distinct += bin_distinct;
            }
        }
    }
    return distinct;
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
    if (0) {
        FILE *fp;
        char kmer_file[100];
        sprintf(kmer_file,"%s/%zi-mers.%zi.%u.raw",outdir,k,bin,batches[bin]);
        fp = fopen(kmer_file, "wb");
        if (nwords == 1)
            fwrite(kmerBuf1[bin].data(), kmerSize, kmerBuf1[bin].size(), fp);
        else
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
}

void Kmerizer::save() {
    state = SAVE;
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
    if (nwords == 1) uniqify1(bin);
    else uniqify(bin);
}

void Kmerizer::mergeBin(size_t bin) {
    // create new bitmap indexes for the kmers (binary encoded) and their frequencies (range encoded)
    // open each batch of files (mmap)
    // review notes on hybrid merge
}