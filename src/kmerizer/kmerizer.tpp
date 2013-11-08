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

    uniformBins = true;
    uint32_t capacity = 100000;
    this->reservedSpace = maximem / kmerSize / NBINS;
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
    float ratio = (float)vacancies/(float)entries;
    // fprintf(stderr,"resizeBins1() entries: %u vacancies: %u reservedSpace: %u\n",entries,vacancies,reservedSpace);

    for(size_t bin=0;bin<NBINS;bin++) {
        kmerBufSize[bin] = kmerBuf1[bin].size() + 1 + floor(kmerBuf1[bin].size()*ratio);
        // fprintf(stderr,"kmerBufSize[%zi]=%u\n",bin,kmerBufSize[bin]);
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
    if (state == SAVE || kmerBuf1[bin].size() > 0.9*(float)kmerBufSize[bin]) {
        vector<uint32_t> tally;
        // output distinct kmers
        // if (kmerBufTally[bin] > 0) {
        if (kmerBuf1[bin].size() > 0) {
            binSortCount1(kmerBuf1[bin],tally);
            // sortCount1(kmerBuf1[bin],tally);
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
        vector<vector<kword_t> > keys;// [hbins];
        keys.resize(hbins);
        
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
    if (1) {
        char kmer_file[100];
        sprintf(kmer_file,"%s/%zi-mers.%zi.%u.bsi",outdir,k,bin,batches[bin]);
        BitSlicedIndex<kword_t> *bsi = new BitSlicedIndex<kword_t>(nwords);
        if (nwords == 1) {
            kword_t *buf = kmerBuf1[bin].data();
            for(int i = 0; i < kmerBuf1[bin].size(); i++)
                bsi->append(buf + i);
        }
        else
            for(size_t i=0;i<kmerBufTally[bin];i++)
                bsi->append(kmerBuf[bin] + i*nwords);
        bsi->saveIndex(kmer_file);
        delete bsi;

        char kmerCount_file[100];
        sprintf(kmerCount_file,"%s/%zi-mers.%zi.%u.rei",outdir,k,bin,batches[bin]);
        RangeEncodedIndex *rei = new RangeEncodedIndex(tally);
        rei->saveIndex(kmerCount_file);
        delete rei;
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
        // mergeBin(bin);
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
    if (batches[bin] == 0) return;
    // create new bitmap indexes for the kmers (binary encoded) and their frequencies (range encoded)
    // review notes on hybrid merge
    if (batches[bin]==1) {
        // nothing to merge - only one batch for this bin
        // move the .bsi and .rei files from .1.rei to .rei
        char old_fname[100], new_fname[100];
        sprintf(old_fname,"%s/%zi-mers.%zi.%u.bsi",outdir,k,bin,batches[bin]);
        sprintf(new_fname,"%s/%zi-mers.%zi.bsi",outdir,k,bin);
        int result = rename(old_fname, new_fname);

        sprintf(old_fname,"%s/%zi-mers.%zi.%u.rei",outdir,k,bin,batches[bin]);
        sprintf(new_fname,"%s/%zi-mers.%zi.rei",outdir,k,bin);
        result = rename(old_fname, new_fname);

        return;
    }
    // need to merge the batches
    // load each from disk
    BitSlicedIndex<kword_t> **bsi;
    RangeEncodedIndex **rei;
    size_t *batchOffset;
    size_t nBatches = batches[bin];

    bsi = (BitSlicedIndex<kword_t> **) malloc(nBatches*sizeof(BitSlicedIndex<kword_t>*));
    rei = (RangeEncodedIndex **) malloc(nBatches*sizeof(RangeEncodedIndex*));
    batchOffset = (size_t *) calloc(nBatches,sizeof(size_t));
    for(size_t i=0;i<nBatches;i++) {
        char batchFname[100];
        sprintf(batchFname,"%s/%zi-mers.%zi.%zi.bsi",outdir,k,bin,i+1);
        bsi[i] = new BitSlicedIndex<kword_t>(batchFname);
        sprintf(batchFname,"%s/%zi-mers.%zi.%zi.rei",outdir,k,bin,i+1);
        rei[i] = new RangeEncodedIndex(batchFname);
    }
    // for the merged bit sliced index
    BitSlicedIndex<kword_t> *mergedBsi = new BitSlicedIndex<kword_t>(nwords);
    vector<uint32_t> mergedCounts;
    mergedCounts.reserve(2 * bsi[0]->size());
    fprintf(stderr,"mergedCounts.capacity() %zi\n",mergedCounts.capacity());
    // diy priority_queue to manage the merge
    kword_t *kmerPairs;
    kmerPairs = (kword_t*) calloc(nBatches*(nwords+1), sizeof(kword_t)); // the extra word holds the batch number
    // load a kmerPair from each batch
    size_t offset=0;
    for(size_t i=0;i<nBatches;i++) {
        if (bsi[i]->decode(0,kmerPairs + offset)) {
            offset += nwords;
            kmerPairs[offset]=(kword_t)i;
            offset++;
        }
        else {
            fprintf(stderr,"batch %zi is empty!\n",i);
            exit(3);
        }
    }
    offset=0;
    
    kword_t *maxKmer, *kmer;
    int nbits = 12; // tune me
    kword_t lutSize = 1 << nbits;
    uint64_t *inSet; // diy bitvector to mark set bits
    uint32_t *counts; // direct hash of merged counts
    uint32_t v;
    inSet = (uint64_t*) calloc(lutSize >> 6, sizeof(uint64_t));
    counts = (uint32_t*) calloc(lutSize, sizeof(uint32_t));
    maxKmer = (kword_t*) malloc(nwords*sizeof(kword_t));
    kmer = (kword_t*) malloc(nwords*sizeof(kword_t));
    while (nBatches > 1) {
        // for(size_t i=0;i<nBatches;i++) {
        //     fprintf(stderr,"unsorted %zi ",i);
        //     printBits(sizeof(kword_t), kmerPairs + i*(nwords+1));
        //     fprintf(stderr," (%llu) batch: %llu\n",kmerPairs[i*(nwords+1)],kmerPairs[1+i*(nwords+1)]);
        // }
        // sort the kmerPairs
        if      (nwords == 1) sort((kmer2_t*)kmerPairs, (kmer2_t*)(kmerPairs + (nwords+1)*nBatches));
        else if (nwords == 2) sort((kmer3_t*)kmerPairs, (kmer3_t*)(kmerPairs + (nwords+1)*nBatches));
        else if (nwords == 3) sort((kmer4_t*)kmerPairs, (kmer4_t*)(kmerPairs + (nwords+1)*nBatches));
        else if (nwords == 4) sort((kmer5_t*)kmerPairs, (kmer5_t*)(kmerPairs + (nwords+1)*nBatches));
        else if (nwords == 5) sort((kmer6_t*)kmerPairs, (kmer6_t*)(kmerPairs + (nwords+1)*nBatches));
        else if (nwords == 6) sort((kmer7_t*)kmerPairs, (kmer7_t*)(kmerPairs + (nwords+1)*nBatches));
        else if (nwords == 7) sort((kmer8_t*)kmerPairs, (kmer8_t*)(kmerPairs + (nwords+1)*nBatches));
        else if (nwords == 8) sort((kmer9_t*)kmerPairs, (kmer9_t*)(kmerPairs + (nwords+1)*nBatches));

        // fprintf(stderr,"nBatches: %zi\n",nBatches);
        // for(int i=0;i<nBatches;i++) {
        //     for(int k=0;k<nwords;k++) {
        //         fprintf(stderr,"i: %i, batch: %llu, kword[%i]: ",
        //         i,kmerPairs[i*(nwords+1)+nwords],k);
        //         printBits(sizeof(kword_t), kmerPairs + i*(nwords+1)+k);
        //         fprintf(stderr,"\n");
        //     }
        // }
        // exit(0);
        // take the minimum value (pointer to beginning of list)
        kword_t *minKmer = kmerPairs;
        // maximum value that will fit in the lookup table + 1
        for(int i=0;i<nwords;i++) maxKmer[i]=minKmer[i];
        if ((kword_t)~(kword_t)0 - minKmer[nwords-1] < lutSize)
            maxKmer[nwords-1] = (kword_t)~(kword_t)0;
        else
            maxKmer[nwords-1] += lutSize;

        vector<kword_t> remainingBatches;
        for(int i=0;i<nBatches;i++) {
            kword_t batch = kmerPairs[i*(nwords+1) + nwords];
            bool check = bsi[batch]->decode(batchOffset[batch],kmer);
            while (check) {
                // for(int j=0;j<nwords;j++) fprintf(stderr,"%i %i kmer:%llu min:%llu max:%llu\n",i,j,kmer[j],minKmer[j],maxKmer[j]);
                if (kmerless(kmer,maxKmer,nwords)) {
                    if (kmerless(kmer,minKmer,nwords)) {
                        fprintf(stderr,"kmers out of order\nkmer:");
                        printBits(sizeof(kword_t), kmer);
                        fprintf(stderr,"\n min:");
                        printBits(sizeof(kword_t), minKmer);
                        fprintf(stderr,"\n");
                    }
                    kword_t koffset = kmer[nwords-1] - minKmer[nwords-1];
                    // fprintf(stderr,"i %i batch %llu koffset %llu\n",i,batch,koffset);
                    inSet[koffset >> 6] |= (kword_t)1 << (koffset & 63);
                    if (rei[batch]->decode(batchOffset[batch], &v))
                        counts[koffset] += v;
                    else
                        { fputs("error decoding from rei",stderr); exit(3); }
                    batchOffset[batch]++;
                    check = bsi[batch]->decode(batchOffset[batch],kmer);
                }
                else {
                    // fprintf(stderr,"kmer ge max\n");
                    check = false;
                    remainingBatches.push_back(batch);
                }
            }
        }
        if (nBatches > remainingBatches.size())
            fprintf(stderr,"exhausted a batch nBatches now %zi\n",nBatches);
        nBatches = remainingBatches.size();

        // populate the kmerPairs array for each of the remaining batches
        for(int i=0;i<nBatches;i++) {
            kword_t batch = remainingBatches[i];
            bsi[batch]->decode(batchOffset[batch],kmerPairs + i*(nwords+1));
            kmerPairs[i*(nwords+1) + nwords] = batch;
        }

        // loop over the set bits and reconstruct (merged) kmers
        // encode the kmers and save the counts for after the loop in mergedCounts.push_back(counts[offset])
        // because the range encoded bitmap index needs all the data before it can do anything
        int iter = lutSize >> 6;
        for(int i=0;i<iter;i++) {
            uint64_t bits = inSet[i];
            inSet[i]=(uint64_t)0;
            int oi = (i<<6) - 1; // - 1 because ffs(bits) is 1 based
            while (bits) {
                offset = ffs(bits) + oi;
                mergedCounts.push_back(counts[offset]);
                counts[offset] = 0; // reset to 0 for next time
                maxKmer[nwords-1] = minKmer[nwords-1] + offset;
                mergedBsi->append(maxKmer);
                bits &= bits-1;
            }
        }
    }
    if (nBatches==1) {
        kword_t batch = kmerPairs[nwords];
        while (bsi[batch]->decode(batchOffset[batch],kmer)) {
            mergedBsi->append(kmer);
            if (rei[batch]->decode(batchOffset[batch], &v))
                mergedCounts.push_back(v);
            batchOffset[batch]++;
        }
    }
    fprintf(stderr,"after merge mergedCounts.size() %zi mergedBsi.size() %zi\n",mergedCounts.size(),mergedBsi->size());
    
    char merged_fname[100];
    sprintf(merged_fname,"%s/%zi-mers.%zi.bsi",outdir,k,bin);
    fprintf(stderr,"mergedBsi->saveIndex(%s)\n",merged_fname);
    mergedBsi->saveIndex(merged_fname);
    delete mergedBsi;
    sprintf(merged_fname,"%s/%zi-mers.%zi.rei",outdir,k,bin);
    fprintf(stderr,"RangeEncodedIndex(mergedCounts) size=%zi\n",mergedCounts.size());
    RangeEncodedIndex *mergedRei = new RangeEncodedIndex(mergedCounts);
    fprintf(stderr,"mergedRei->saveIndex(%s) mergedRei->size() %zi\n",merged_fname,mergedRei->size());
    mergedRei->saveIndex(merged_fname);
    delete mergedRei;

    for(int i=0;i<batches[bin];i++) {
        delete bsi[i];
        delete rei[i];
        char batchFname[100];
        sprintf(batchFname,"%s/%zi-mers.%zi.%i.bsi",outdir,k,bin,i+1);
        if (remove(batchFname)) { fprintf(stderr, "error removing %s\n",batchFname); exit(3); }
        sprintf(batchFname,"%s/%zi-mers.%zi.%i.rei",outdir,k,bin,i+1);
        if (remove(batchFname)) { fprintf(stderr, "error removing %s\n",batchFname); exit(3); }
    }
}
