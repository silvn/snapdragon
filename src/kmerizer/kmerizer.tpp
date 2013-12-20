// #include <bitset>
#include <algorithm> // stl sort
// #include <cstdlib> // qsort()
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

    this->kmask       = (kword_t)~(kword_t)0;
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

    memset(filtered,0,NBINS);
}

int Kmerizer::allocate(const size_t maximem) {
    fprintf(stderr, "Kmerizer::allocate(%zi)", maximem);
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);

    memset(batches, 0, sizeof(uint32_t) * NBINS);

    if (nwords>1)
        memset(kmerBufTally, 0, sizeof(uint32_t) * NBINS);

    this->reservedSpace = maximem / kmerSize / NBINS;
    uint32_t capacity = reservedSpace; //100000;

    if (nwords == 1)
        uniformBins = true;
    else
        capacity = reservedSpace;
     
    if (reservedSpace > capacity)
        reservedSpace -= capacity;
    else
        reservedSpace=0;
    for (size_t i = 0; i < NBINS; i++) {
        kmerBufSize[i] = capacity;
        if (nwords == 1)
            kmerBuf1[i].reserve(capacity);
        else
            kmerBuf[i] = (kword_t *) calloc(capacity, kmerSize);
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
    while (a.a[i] == b.a[i]) {
        i++;
        if (i==s) return false; // ==
    }
    return a.a[i] < b.a[i];
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
bool operator<(const kmer9_t& a, const kmer9_t& b) {
    int i=0, s=9;
    while (*(&a.a+i) == *(&b.a+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(&a.a+i) < *(&b.a+i);
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
        rcpack[nwords-2] |= *bin << (56+lshift);
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

void Kmerizer::unpack(kword_t *kmer, size_t bin, char *kmerStr) {
    // space already allocated for kmerStr of length k+4
    static const char lut[4] = {'A','C','G','T'}; // check this
    kmerStr[0] = lut[bin >> 6];
    kmerStr[1] = lut[(bin >> 4) & 3];
    kmerStr[2] = lut[(bin >> 2) & 3];
    kmerStr[3] = lut[bin & 3];
    if (nwords == 1)
        for(int i=0;i<k;i++)
            kmerStr[4+i] = lut[(kmer[0] >> 2*(k-i-1)) & 3];
    else {
        for(int w=0;w<nwords-1;w++)
            for(int i=0;i<k;i++)
                kmerStr[4 + i + w*32] = lut[(kmer[w] >> 2*(31-i)) & 3];
        // last word
        int kmod = k%32;
        for(int i=0;i<kmod;i++)
            kmerStr[4 + i + (nwords-1)*32] = lut[(kmer[nwords-1] >> 2*(kmod-i-1)) & 3];
    }
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
            if (uniformBins) {
                resizeBins1();
                uniformBins = false;
            }
            else {
                for(size_t i=0;i<NBINS;i++)
                    tp.schedule( boost::bind( &Kmerizer::uniqify1, this, i ) );
                tp.wait();
            }
        }

        for (int i=k+4; i<length;i++) {
            // bin = nextKmer1(&packed, bin, seq[i]);
            bin = ((bin << 2) | (packed >> shiftlastby)) & 255;
            packed = ((packed << 2) | lut[(seq[i] >> 1) & 3]) & kmask;
            kmer = canonicalize1(&packed,&rcpack,&bin);
            // insertKmer1(bin, kmer);
            kmerBuf1[bin].push_back(*kmer);
            if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
                if (uniformBins) {
                    resizeBins1();
                    uniformBins = false;
                }
                else {
                    for(size_t b=0;b<NBINS;b++)
                        tp.schedule( boost::bind( &Kmerizer::uniqify1, this, b ) );
                    tp.wait();
                }
            }
        }
    }
    else {
        // insertKmer1(bin, &packed);
        kmerBuf1[bin].push_back(packed);
        if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
            if (uniformBins) {
                resizeBins1();
                uniformBins = false;
            }
            else {
                for(size_t b=0;b<NBINS;b++)
                    tp.schedule( boost::bind( &Kmerizer::uniqify1, this, b ) );
                tp.wait();
            }
        }
        for (int i=k+4; i<length;i++) {
            // bin = nextKmer1(&packed, bin, seq[i]);
            bin = ((bin << 2) | (packed >> shiftlastby)) & 255;
            packed = ((packed << 2) | lut[(seq[i] >> 1) & 3]) & kmask;
            // insertKmer1(bin, &packed);
            kmerBuf1[bin].push_back(packed);
            if (kmerBuf1[bin].size() == kmerBufSize[bin]) {
                if (uniformBins) {
                    resizeBins1();
                    uniformBins = false;
                }
                else {
                    for(size_t b=0;b<NBINS;b++)
                        tp.schedule( boost::bind( &Kmerizer::uniqify1, this, b ) );
                    tp.wait();
                }
            }
        }
    }
}

/*
    count          x
    -----   =  ---------
    entries    vacancies

e*x = c*v
x = c*v/e
    
*/
void Kmerizer::resizeBins1() {
    // measure the total amount of free space
    uint32_t entries=0;
    uint32_t vacancies=reservedSpace*NBINS;
    for(size_t bin=0;bin<NBINS;bin++) {
        entries += kmerBuf1[bin].size();
        vacancies += kmerBufSize[bin] - kmerBuf1[bin].size();
    }
    double ratio = (double)vacancies/(double)entries;

    for(size_t bin=0;bin<NBINS;bin++) {
        kmerBufSize[bin] = kmerBuf1[bin].size() + 1 + floor(ratio*(double)kmerBuf1[bin].size());
        if (kmerBufSize[bin] < kmerBuf1[bin].capacity()) {
            // make a copy and swap
            vector<kword_t> v (kmerBuf1[bin]);
            kmerBuf1[bin].swap(v);
        }
        else
            kmerBuf1[bin].reserve(kmerBufSize[bin]);
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

bool kmerless(kword_t *a, kword_t *b, int s) {
    int i=0;
    while (*(a+i) == *(b+i)) {
        i++;
        if (i==s) return false; // ==
    }
    return *(a+i) < *(b+i);
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

    // if sampling phase then create an ordered map to count frequent k-mers
    // copy the first sampleSize kmers into another vector
    // do the binSortCount1()
    // identify the most frequent LUTSize kmers
    // copy them into the LUT

    // if (state == SAVE || kmerBufTally[bin] > 0.9*(float)kmerBufSize[bin]) {
    if (state == SAVE || kmerBuf1[bin].size() > floor(0.9*(float)kmerBufSize[bin])) {
        vector<uint32_t> tally;
        // output distinct kmers
        // if (kmerBufTally[bin] > 0) {
        if (kmerBuf1[bin].size() > 0) {
            // fprintf(stderr,"%zi\t%zi\n",bin,kmerBuf1[bin].size());
            // binSortCount1(kmerBuf1[bin],tally);
            // sortCountKmers(kmerBuf1[bin],tally);
            // sortCount1(kmerBuf1[bin],tally);
            batches[bin]++;

            // writeBatch(bin, tally);
            sortwriteBatch(bin);
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

void Kmerizer::sortCount1(vector<kword_t> &kb, BitSlicedIndex<kword_t> *kmerIdx, LCBitSlicedIndex<uint32_t> *countIdx) {
    sort(kb.begin(), kb.end());
    // uniq -c
    kword_t kmer = kb[0];
    uint32_t tally = 1; // first kmer
    for(vector<kword_t>::iterator it = kb.begin()+1; it < kb.end(); it++) {
        if(*it == kmer)
            tally++;
        else {
            kmerIdx->append(&kmer);
            countIdx->append(&tally);
            tally=1;
            kmer = *it;
        }
    }
    kmerIdx->append(&kmer);
    countIdx->append(&tally);
}

void Kmerizer::binSortCount1(vector<kword_t> &kb, BitSlicedIndex<kword_t> *kmerIdx, LCBitSlicedIndex<uint32_t> *countIdx) {
    if (kb.size() < (1 << 18)) {
        sortCount1(kb,kmerIdx,countIdx);
        return;
    }
    int hbits = 12;
    int hbins = 1 << hbits; // 4096
    uint32_t *histogram = (uint32_t *) calloc(hbins, sizeof(uint32_t));

    int shiftby = 2*k - hbits;
    if (shiftby <= 0) {
        // direct hashing
        for(size_t i=0;i<kb.size();i++)
            histogram[kb[i]]++;

        for(uint32_t b = 0; b < hbins; b++)
            if (histogram[b] > 0) {
                kword_t kmerized_b = (kword_t)b;
                kmerIdx->append(&kmerized_b);
                countIdx->append(histogram + b);
            }
    }
    else {
        vector<vector<kword_t> > keys;// [hbins];
        keys.resize(hbins);
        // count the number of kmers in each bin
        for(vector<kword_t>::iterator it = kb.begin(); it < kb.end(); it++)
            histogram[*it >> shiftby]++;
        
        // allocate space for each bin
        for(int i=0; i<hbins; i++)
            if (histogram[i] > 0)
                keys[i].reserve(histogram[i]);

        // loop over the kmers again and copy them into their bins
        for(vector<kword_t>::iterator it = kb.begin(); it < kb.end(); it++)
            keys[*it >> shiftby].push_back(*it);
        
        for(uint32_t b = 0; b < hbins; b++)
            if (keys[b].size() > 0)
                sortCount1(keys[b],kmerIdx,countIdx);
    }
    free(histogram);
}

void Kmerizer::binSortCount1(vector<kword_t> &kb, vector<uint32_t> &tally) {
    if (kb.size() < (1 << 18)) {
        sortCount1(kb,tally);
        return;
    }
    int hbits = 12;
    int hbins = 1 << hbits; // 4096
    uint32_t *histogram = (uint32_t *) calloc(hbins, sizeof(uint32_t));

    int shiftby = 2*k - hbits;
    if (shiftby <= 0) {
        // direct hashing
        for(size_t i=0;i<kb.size();i++)
            histogram[kb[i]]++;

        kb.clear();
        for(uint32_t b = 0; b < hbins; b++)
            if (histogram[b] > 0) {
                kb.push_back((kword_t) b);
                tally.push_back(histogram[b]);
            }
    }
    else {
        vector<vector<kword_t> > keys;// [hbins];
        keys.resize(hbins);
        // count the number of kmers in each bin
        // for(size_t i=0;i<kb.size();i++)
        for(vector<kword_t>::iterator it = kb.begin(); it < kb.end(); it++)
            histogram[*it >> shiftby]++;
        
        // allocate space for each bin
        for(int i=0; i<hbins; i++)
            if (histogram[i] > 0)
                keys[i].reserve(histogram[i]);

        // loop over the kmers again and copy them into their bins
        // for(uint32_t i=0;i<kb.size();i++)
        for(vector<kword_t>::iterator it = kb.begin(); it < kb.end(); it++)
            keys[*it >> shiftby].push_back(*it);
        
        kb.clear();
        for(uint32_t b = 0; b < hbins; b++) {
            if (keys[b].size() > 0) {
                sortCount1(keys[b],tally);
                kb.insert(kb.end(),keys[b].begin(),keys[b].end());
            }
        }
    }
    free(histogram);
}


void Kmerizer::writeBatch(size_t bin, vector<uint32_t> &tally) {
    // if (bin != 0) return;
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
        fwrite(tally.data(), sizeof(uint32_t), tally.size(), fp);
        fclose(fp);
    }
    if (0) {
        FILE *fp;
        char kmer_file[100];
        sprintf(kmer_file,"%s/%zi-mers.%zi.%u.txt",outdir,k,bin,batches[bin]);
        fp = fopen(kmer_file, "w");
        for(int i=0;i<kmerBuf1[bin].size();i++)
            fprintf(fp,"%llu\t%u\n",kmerBuf1[bin][i], tally[i]);
        fclose(fp);
    }
    if (1) {
        char kmer_file[100];
        sprintf(kmer_file,"%s/%zi-mer.seq.%zi.%u",outdir,k,bin,batches[bin]);
        BitSlicedIndex<kword_t> *bsi = new BitSlicedIndex<kword_t>(nwords);
        char kmerCount_file[100];
        sprintf(kmerCount_file,"%s/%zi-mer.count.%zi.%u",outdir,k,bin,batches[bin]);
        LCBitSlicedIndex<uint32_t> *lcbsi = new LCBitSlicedIndex<uint32_t>();
        uint32_t *cbuf = tally.data();
        if (nwords == 1) {
            kword_t *buf = kmerBuf1[bin].data();
            size_t n_vals = kmerBuf1[bin].size();
            for(int i = 0; i < n_vals; i++) {
                bsi->append(buf + i);
                lcbsi->append(cbuf + i);
            }
        }
        else
            for(size_t i=0;i<kmerBufTally[bin];i++) {
                bsi->append(kmerBuf[bin] + i*nwords);
                lcbsi->append(cbuf + i);
            }
        // bsi->saveIndex(kmer_file);
        // lcbsi->saveIndex(kmerCount_file);
        delete bsi;
        delete lcbsi;

        // RangeEncodedIndex *rei = new RangeEncodedIndex(tally);
        // rei->saveIndex(kmerCount_file);
        // delete rei;
    }
}
void Kmerizer::sortwriteBatch(size_t bin) {
    // if (bin != 0) return;
    if (1) {
        char kmer_file[100];
        sprintf(kmer_file,"%s/%zi-mer.seq.%zi.%u",outdir,k,bin,batches[bin]);
        BitSlicedIndex<kword_t> *bsi = new BitSlicedIndex<kword_t>(nwords);
        char kmerCount_file[100];
        sprintf(kmerCount_file,"%s/%zi-mer.count.%zi.%u",outdir,k,bin,batches[bin]);
        LCBitSlicedIndex<uint32_t> *lcbsi = new LCBitSlicedIndex<uint32_t>();

        binSortCount1(kmerBuf1[bin],bsi,lcbsi);

        bsi->saveIndex(kmer_file);
        lcbsi->saveIndex(kmerCount_file);
        delete bsi;
        delete lcbsi;
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

    // clean up memory
    if (nwords == 1)
        for (size_t i = 0; i < NBINS; i++) kmerBuf1[i].clear();
    else
        for (size_t i = 0; i < NBINS; i++) free(kmerBuf[i]);
    
    state = MERGE;
    // merge the batches
    fprintf(stderr,"starting merge\n");
    // mergeBin(0);
    if (nwords == 1)
        for(size_t bin=0; bin < NBINS; bin++)
           tp.schedule( boost::bind( &Kmerizer::mergeBin1, this, bin ) );
    else
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
    // fprintf(stderr,"mergeBin(%zi) batches: %zi\n",bin,batches[bin]);
    if (batches[bin] == 0) return;
    // create new bitmap indexes for the kmers (binary encoded) and their frequencies (range encoded)
    if (batches[bin]==1) {
        // nothing to merge - only one batch for this bin
        // rename the index files
        char old_fname[100], new_fname[100];
        sprintf(old_fname,"%s/%zi-mer.seq.%zi.%u",outdir,k,bin,batches[bin]);
        sprintf(new_fname,"%s/%zi-mer.seq.%zi",outdir,k,bin);
        int result = rename(old_fname, new_fname);

        sprintf(old_fname,"%s/%zi-mer.count.%zi.%u",outdir,k,bin,batches[bin]);
        sprintf(new_fname,"%s/%zi-mer.count.%zi",outdir,k,bin);
        result = rename(old_fname, new_fname);

        return;
    }
    // need to merge the batches
    // load each from disk
    BitSlicedIndex<kword_t> **bsi;
    LCBitSlicedIndex<uint32_t> **lcbsi;
    size_t *batchOffset;
    size_t nBatches = batches[bin];

    bsi = (BitSlicedIndex<kword_t> **) malloc(nBatches*sizeof(BitSlicedIndex<kword_t>*));
    lcbsi = (LCBitSlicedIndex<uint32_t> **) malloc(nBatches*sizeof(LCBitSlicedIndex<uint32_t>*));
    batchOffset = (size_t *) calloc(nBatches,sizeof(size_t));
    for(size_t i=0;i<nBatches;i++) {
        char batchFname[100];
        sprintf(batchFname,"%s/%zi-mer.seq.%zi.%zi",outdir,k,bin,i+1);
        bsi[i] = new BitSlicedIndex<kword_t>(batchFname);
        sprintf(batchFname,"%s/%zi-mer.count.%zi.%zi",outdir,k,bin,i+1);
        lcbsi[i] = new LCBitSlicedIndex<uint32_t>(batchFname);
    }
    // for the merged bit sliced index
    BitSlicedIndex<kword_t> *mergedBsi = new BitSlicedIndex<kword_t>(nwords);
    LCBitSlicedIndex<uint32_t> *mergedlcBsi = new LCBitSlicedIndex<uint32_t>();

    kword_t *batchKmers;
    batchKmers = (kword_t*) calloc(nBatches*nwords, sizeof(kword_t));

    // load a kmerPair from each batch
    for(size_t i=0;i<nBatches;i++) {
        if (! bsi[i]->decode(0,batchKmers + i*nwords)) {
            fprintf(stderr,"batch %zi is empty!\n",i);
            exit(3);
        }
    }
    
    kword_t *minKmer, *maxKmer, *kmer;
    int nbits = 12; // tune me
    kword_t lutSize = 1 << nbits;
    uint64_t *inSet; // diy bitvector to mark set bits
    uint32_t *counts; // direct hash of merged counts
    uint32_t v;
    inSet = (uint64_t*) calloc(lutSize >> 6, sizeof(uint64_t));
    counts = (uint32_t*) calloc(lutSize, sizeof(uint32_t));
    kmer = (kword_t*) malloc(nwords*sizeof(kword_t));
    maxKmer = (kword_t*) malloc(nwords*sizeof(kword_t));
    vector<size_t> remainingBatches, todo;
    for(size_t b=0;b<nBatches;b++)
        remainingBatches.push_back(b);
    
    while (remainingBatches.size() > 1) {
        // select minKmer from the batchKmers
        vector<size_t>::iterator bi = remainingBatches.begin();
        minKmer = batchKmers + *bi*nwords;
        bi++;
        while (bi != remainingBatches.end()) {
            if(kmerless(batchKmers + *bi*nwords, minKmer, nwords))
                minKmer = batchKmers + *bi*nwords;
            bi++;
        }
        
        // maximum value that will fit in the lookup table + 1
        for(int i=0;i<nwords;i++) maxKmer[i]=minKmer[i];
        if ((kword_t)~(kword_t)0 - minKmer[nwords-1] < lutSize)
            maxKmer[nwords-1] = (kword_t)~(kword_t)0;
        else
            maxKmer[nwords-1] += lutSize;

        for(bi = remainingBatches.begin(); bi != remainingBatches.end(); bi++) {
            bool check = bsi[*bi]->decode(batchOffset[*bi],kmer);
            while (check) {
                if (kmerless(kmer,maxKmer,nwords)) {
                    kword_t koffset = kmer[nwords-1] - minKmer[nwords-1];
                    inSet[koffset >> 6] |= ((kword_t)1 << (koffset & 63));
                    // if (rei[*bi]->decode(batchOffset[*bi], &v))
                    // fprintf(stderr,"cbsi[%zi]->decode(%zi,v)\n",*bi,batchOffset[*bi]);
                    if (lcbsi[*bi]->decode(batchOffset[*bi], &v))
                        counts[koffset] += v;
                    else
                        { fputs("error decoding from lcbsi",stderr); exit(3); }
                    batchOffset[*bi]++;
                    check = bsi[*bi]->decode(batchOffset[*bi],kmer);
                }
                else {
                    check = false;
                    todo.push_back(*bi);
                }
            }
        }
        // loop over the set bits and reconstruct (merged) kmers
        // encode the kmers and save the counts for after the loop in mergedCounts.push_back(counts[offset])
        // because the range encoded bitmap index needs all the data before it can do anything
        // fprintf(stderr,"reconstruct kmers from lut\n");
        int iter = lutSize >> 6;
        for(int i=0;i<iter;i++) {
            uint64_t bits = inSet[i];
            if (bits) {
                inSet[i]=(uint64_t)0;
                int oi = (i<<6) - 1; // - 1 because ffs(bits) is 1 based
                while (bits) {
                    int offset = __builtin_ffsll(bits) + oi;
                    // mergedCounts.push_back(counts[offset]);
                    mergedlcBsi->append(counts + offset);
                    counts[offset] = 0; // reset to 0 for next time
                    maxKmer[nwords-1] = minKmer[nwords-1] + offset;
                    mergedBsi->append(maxKmer);

                    bits &= bits-1;
                }
            }
        }
        
        // finished a batch?
        if (remainingBatches.size() > todo.size())
            swap(remainingBatches,todo);
        todo.clear();
        
        // populate the batchKmers array for each of the remaining batches
        for(bi = remainingBatches.begin(); bi!=remainingBatches.end();bi++)
            bsi[*bi]->decode(batchOffset[*bi],batchKmers + *bi*nwords);
    }
    if (remainingBatches.size()==1) {
        size_t batch = remainingBatches[0];
        size_t bo = batchOffset[batch];
        while (bsi[batch]->decode(bo,kmer)) {
            mergedBsi->append(kmer);
            if (lcbsi[batch]->decode(bo, &v))
                mergedlcBsi->append(&v);
            // if (rei[batch]->decode(bo, &v))
            //     mergedCounts.push_back(v);
            bo++;
        }
    }
    for(int i=0;i<batches[bin];i++) {
        delete bsi[i];
        delete lcbsi[i];
        char batchFname[100];
        sprintf(batchFname,"%s/%zi-mer.seq.%zi.%i",outdir,k,bin,i+1);
        if (remove(batchFname)) { fprintf(stderr, "error removing %s\n",batchFname); exit(3); }
        sprintf(batchFname,"%s/%zi-mer.count.%zi.%i",outdir,k,bin,i+1);
        if (remove(batchFname)) { fprintf(stderr, "error removing %s\n",batchFname); exit(3); }
    }
    
    char merged_fname[100];
    sprintf(merged_fname,"%s/%zi-mer.seq.%zi",outdir,k,bin);
    mergedBsi->saveIndex(merged_fname);
    delete mergedBsi;

    sprintf(merged_fname,"%s/%zi-mer.count.%zi",outdir,k,bin);
    mergedlcBsi->saveIndex(merged_fname);
    delete mergedlcBsi;
}

struct kcb_type {
    kword_t kmer;
    uint32_t tally;
    size_t batch;
};
class kcb_cmp {
public:
    kcb_cmp() { }
    bool operator() (const kcb_type& lhs, const kcb_type& rhs) const {
        return (lhs.kmer > rhs.kmer);
    }
};
typedef priority_queue<kcb_type, vector<kcb_type>, kcb_cmp> kcb_pq_type;

void Kmerizer::mergeBin1(size_t bin) {
    if (batches[bin] == 0) return;
    // create new bitmap indexes for the kmers and their frequencies
    if (batches[bin]==1) {
        // nothing to merge - only one batch for this bin
        // rename the index files
        char old_fname[100], new_fname[100];
        sprintf(old_fname,"%s/%zi-mer.seq.%zi.%u",outdir,k,bin,batches[bin]);
        sprintf(new_fname,"%s/%zi-mer.seq.%zi",outdir,k,bin);
        int result = rename(old_fname, new_fname);

        sprintf(old_fname,"%s/%zi-mer.count.%zi.%u",outdir,k,bin,batches[bin]);
        sprintf(new_fname,"%s/%zi-mer.count.%zi",outdir,k,bin);
        result = rename(old_fname, new_fname);

        return;
    }
    // need to merge the batches
    // load each from disk
    BitSlicedIndex<kword_t> **bsi;
    LCBitSlicedIndex<uint32_t> **lcbsi;

    size_t *batchOffset;
    size_t nBatches = batches[bin];

    bsi = (BitSlicedIndex<kword_t> **) malloc(nBatches*sizeof(BitSlicedIndex<kword_t>*));
    lcbsi = (LCBitSlicedIndex<uint32_t> **) malloc(nBatches*sizeof(LCBitSlicedIndex<uint32_t>*));

    batchOffset = (size_t *) calloc(nBatches,sizeof(size_t));
    for(size_t i=0;i<nBatches;i++) {
        char batchFname[100];
        sprintf(batchFname,"%s/%zi-mer.seq.%zi.%zi",outdir,k,bin,i+1);
        bsi[i] = new BitSlicedIndex<kword_t>(batchFname);
        sprintf(batchFname,"%s/%zi-count.seq.%zi.%zi",outdir,k,bin,i+1);
        lcbsi[i] = new LCBitSlicedIndex<uint32_t>(batchFname);
    }
    // for the merged bit sliced index
    BitSlicedIndex<kword_t> *mergedBsi = new BitSlicedIndex<kword_t>(1);
    LCBitSlicedIndex<uint32_t> *mergedlcBsi = new LCBitSlicedIndex<uint32_t>();

    // put the first kmer and its frequency from each batch into a priority queue
    kcb_pq_type pq;

    kcb_type kcb_triple;
    kword_t kmer;
    uint32_t tally;
    size_t batch;
    for(size_t i=0;i<nBatches;i++) {
        if (! bsi[i]->decode(0,&(kcb_triple.kmer))) {
            fprintf(stderr,"batch %zi is empty!\n",i);
            exit(3);
        }
        if (! lcbsi[i]->decode(0,&(kcb_triple.tally))) {
            fprintf(stderr,"count batch %zi is empty\n",i);
            exit(3);
        }
        kcb_triple.batch = i;
        pq.push(kcb_triple);
    }

    // iterate until the priority queue is empty
    kcb_triple = pq.top();
    kmer = kcb_triple.kmer;
    tally = 0;
    while (!pq.empty()) {
        kcb_triple = pq.top();
        pq.pop();
        if (kcb_triple.kmer == kmer)
            tally += kcb_triple.tally;
        else {
            mergedBsi->append(&kmer);
            mergedlcBsi->append(&tally);
            kmer = kcb_triple.kmer;
            tally = kcb_triple.tally;
        }
        batch = kcb_triple.batch;
        batchOffset[batch]++;
        if (bsi[batch]->decode(batchOffset[batch],&(kcb_triple.kmer)) &&
            lcbsi[batch]->decode(batchOffset[batch],&(kcb_triple.tally)))
            pq.push(kcb_triple);
    }
    mergedBsi->append(&kmer);
    mergedlcBsi->append(&tally);


    for(int i=0;i<batches[bin];i++) {
        delete bsi[i];
        delete lcbsi[i];
        char batchFname[100];
        sprintf(batchFname,"%s/%zi-mer.seq.%zi.%i",outdir,k,bin,i+1);
        if (remove(batchFname)) { fprintf(stderr, "error removing %s\n",batchFname); exit(3); }
        sprintf(batchFname,"%s/%zi-mer.count.%zi.%i",outdir,k,bin,i+1);
        if (remove(batchFname)) { fprintf(stderr, "error removing %s\n",batchFname); exit(3); }
    }
    
    char merged_fname[100];
    sprintf(merged_fname,"%s/%zi-mer.seq.%zi",outdir,k,bin);
    mergedBsi->saveIndex(merged_fname);
    delete mergedBsi;

    sprintf(merged_fname,"%s/%zi-mer.count.%zi",outdir,k,bin);
    mergedlcBsi->saveIndex(merged_fname);
    delete mergedlcBsi;

}

// unique, distinct, total, max_count from colname (usu. count)
void Kmerizer::stats(char* colname) {
    uint32_t unique=0,distinct=0,total=0,max_count=0;
    for(size_t bin=0; bin < NBINS; bin++) {
        char fname[100];
        sprintf(fname,"%s/%zi-mer.%s.%zi",outdir,k,colname,bin);
        LCBitSlicedIndex<uint32_t> *lcbsi;
        lcbsi = new LCBitSlicedIndex<uint32_t>(fname,true);
        unique += lcbsi->unique();
        distinct += lcbsi->size();
        total += lcbsi->totalCount();
        if (max_count < lcbsi->maxCount()) max_count = lcbsi->maxCount();
    }
    fprintf(stdout,"unique:    %u\ndistinct:  %u\ntotal:     %u\nmax_count: %u\n",unique,distinct,total,max_count);
}

class mycmp2 {
public:
    mycmp2() { }
    bool operator() (const pair <uint32_t, uint32_t>& lhs, const pair <uint32_t, uint32_t>& rhs) const {
        return (lhs.first > rhs.first);
    }
};
typedef priority_queue<pair <uint32_t,uint32_t>, vector<pair <uint32_t,uint32_t> >, mycmp2> pq32_type;

// for each frequency, report the number of kmers
void Kmerizer::histo(char *colname) {
    char fname[100];
    pq32_type pq;
    for(size_t bin=0; bin < NBINS; bin++) {
        sprintf(fname,"%s/%zi-mer.%s.%zi",outdir,k,colname,bin);
        LCBitSlicedIndex<uint32_t> *lcbsi;
        lcbsi = new LCBitSlicedIndex<uint32_t>(fname,true);
        vector<uint32_t>::iterator i1 = lcbsi->dValues.begin();
        vector<uint32_t>::iterator i2 = lcbsi->dFrequency.begin();
        
        while (i1 != lcbsi->dValues.end()) {
            pq.push(make_pair(*i1,*i2));
            i1++;
            i2++;
        }
    }
    uint32_t v = 0;
    uint32_t f = 0;
    while (!pq.empty()) {
        if (v == pq.top().first)
            f += pq.top().second;
        else {
            if (f>0) fprintf(stdout,"%u\t%u\n",v,f);
            v = pq.top().first;
            f = pq.top().second;
        }
        pq.pop();
    }
    if (f>0) fprintf(stdout,"%u\t%u\n",v,f);
}

// results are stored in BitMask
// sequential calls to filter() AND the result set
void Kmerizer::filter(char *colname, uint32_t min, uint32_t max) {
    // fprintf(stderr,"filter(%s, %u, %u)\n",colname, min,max);

    for(size_t bin=0; bin < NBINS; bin++)
        tp.schedule( boost::bind( &Kmerizer::filterBin, this, bin, colname, min, max ) );

    tp.wait();

}

void Kmerizer::filterBin(size_t bin, char* colname, uint32_t min, uint32_t max) {
    char fname[100];
    sprintf(fname,"%s/%zi-mer.%s.%zi",outdir,k,colname,bin);
    // fprintf(stderr,"filterBin(%zi,%s,%u,%u) %s\n",bin,colname,min,max,fname);
    LCBitSlicedIndex<uint32_t> *lcbsi;
    lcbsi = new LCBitSlicedIndex<uint32_t>(fname);
    if (filtered[bin]) // if results exist, AND with this filter
        *BitMask[bin] &= lcbsi->continuousRange(min,max);
    else {
        BitMask[bin] = lcbsi->continuousRange(min,max);
        filtered[bin] = 1;
    }
    // if (BitMask[bin]->getCount() > 0)
    //     fprintf(stderr,"filterBin(%zi,%s,%u,%u) BitMask[bin]->getCount() %zi\n",bin,colname,min,max,BitMask[bin]->getCount());
}

void Kmerizer::dump(vector<char *> &columns) { // write tab delimited text to stdout
    // TODO: named columns (seq, count, <user supplied column names>)
    char fname[100];
    kword_t *kmer;
    char *kmerStr;
    uint32_t count;
    kmer = (kword_t *) malloc(nwords * sizeof(kword_t));
    kmerStr = (char *) malloc(k+5);
    kmerStr[k+4] = '\0';
    for(size_t bin=0;bin < NBINS; bin++) {
        // fprintf(stderr,"dump() bin %zi BitMask[bin]->getCount() %u\n",bin,BitMask[bin]->getCount());
        if (filtered[bin]) {
            if (BitMask[bin]->getCount() > 0) {
                sprintf(fname,"%s/%zi-mer.count.%zi",outdir,k,bin);
                LCBitSlicedIndex<uint32_t> *lcbsi;
                lcbsi = new LCBitSlicedIndex<uint32_t>(fname);
                sprintf(fname,"%s/%zi-mer.seq.%zi",outdir,k,bin);
                BitSlicedIndex<kword_t> *bsi;
                bsi = new BitSlicedIndex<kword_t>(fname);
                // iterate over the set bits and decode from each column
                size_t i=0;
                BitMask[bin]->rewind();
                while (BitMask[bin]->nextSetBit(&i)) {
                    bsi->decode(i,kmer);
                    unpack(kmer,bin,kmerStr);
                    lcbsi->decode(i,&count);
                    fprintf(stdout, "%s\t%u\n",kmerStr,count);
                    i++;
                }
            }
        }
        else {
            // iterate over all positions and decode from each column
            sprintf(fname,"%s/%zi-mer.count.%zi",outdir,k,bin);
            LCBitSlicedIndex<uint32_t> *lcbsi;
            lcbsi = new LCBitSlicedIndex<uint32_t>(fname);
            sprintf(fname,"%s/%zi-mer.seq.%zi",outdir,k,bin);
            BitSlicedIndex<kword_t> *bsi;
            bsi = new BitSlicedIndex<kword_t>(fname);
            for(size_t i=0;i < bsi->size(); i++) {
                bsi->decode(i,kmer);
                unpack(kmer,bin,kmerStr);
                lcbsi->decode(i,&count);
                fprintf(stdout, "%s\t%u\n",kmerStr,count);
            }
        }
    }
}

void Kmerizer::query(const char* seq, const int length) { // extracts kmers, queries indexes, stores results in BitMask
    // the BitMask is initially all 0's and is OR'ed with query results. Subsequent calls only extend the result set
    
}

