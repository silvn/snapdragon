#include "kmerizer.h"
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <bitset>
#include <algorithm> // stl sort
#include <cstdlib> // qsort()
#include <cstdio>  // sprintf()
#include <cstring> // memcpy()
#include <sys/stat.h> // mkdir()
#include <ctime> // clock_t clock() CLOCKS_PER_SEC

Kmerizer::Kmerizer(const size_t k,
                   const size_t threads,
                   const char * outdir,
                   const char   mode)
{
    this->k       = k;
    this->outdir  = new char[strlen(outdir)]; strcpy(this->outdir, outdir);
    this->mode    = mode;
    this->threads = threads;

    this->kmask       = 0xFFFFFFFFFFFFFFFFULL;
    this->shiftlastby = 62;
    this->lshift      = 0;
    this->rshift      = 0;
    if (k % 32 > 0) {
        kmask = (1ULL << (2 * (k % 32))) - 1;
        shiftlastby = 2 * (k % 32) - 2;
        lshift = 2 * (k % 32);
        rshift = 64 - lshift;
    }
    this->nwords     = ((k-1)>>5)+1;
    this->kmerSize   = this->nwords * sizeof(uint32_t);
    this->threadBins = NBINS / this->threads;
    this->state      = READING;
    this->batches    = 0;
}

int Kmerizer::allocate(const size_t maximem) {
    fprintf(stderr, "Kmerizer::allocate(%zi)", maximem);
    clock_t start = clock();
    memset(binTally, 0, sizeof(uint32_t) * NBINS);
    maxKmersPerBin = maximem / kmerSize / NBINS;
    for (size_t i = 0; i < NBINS; i++) {
        kmerBuf[i] = (kword_t *) calloc(maxKmersPerBin, kmerSize);
        if (kmerBuf[i] == NULL) return 1;
    }
    fprintf(stderr," took %f seconds\n",((float)(clock()-start))/CLOCKS_PER_SEC);
    
    return 0;
}

kword_t * Kmerizer::canonicalize(kword_t *packed, kword_t *rcpack) const {
    for (size_t i=0;i<nwords;i++) {
        rcpack[i] = packed[nwords-1-i];
    }
    if (lshift) {
        for (size_t i=0;i<nwords-1;i++) {
            rcpack[i] |= rcpack[i+1] << lshift;
            rcpack[i+1] >>= rshift;
        }
    }
    int cmp=0;
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
    if(cmp > 0)
        return rcpack;
    else
        return packed;
}

void Kmerizer::nextKmer(kword_t* kmer, const char nucl) {
    // shift first kmer by 2 bits
    kmer[0] <<= 2;
    for (size_t w = 1; w < nwords - 1; w++) { // middle (full length) words
        kmer[w-1] |= kmer[w] >> 62; // move the top two bits
        kmer[w] <<= 2; // make room for the next word
    }
    // last word (if not the first)
    if(nwords>1) {
        kmer[nwords-2] |= (kmer[nwords-1] >> shiftlastby);
        kmer[nwords-1] <<= 2;
    }
    kmer[nwords-1] |= twoBit(nucl);
    kmer[nwords-1] &= kmask;
}

void Kmerizer::addSequence(const char* seq, const int length) {
    if (length < k) return;
    kword_t packed[nwords];
    kword_t rcpack[nwords];
    memset(packed, 0, kmerSize);

    for (size_t i = 0; i < k; i++)
        nextKmer(packed,seq[i]);

    kword_t *kmer = packed;
    size_t bin;
    if (mode == CANONICAL)
        kmer = canonicalize(packed,rcpack);
    if (mode == BOTH) {
        kmer = canonicalize(packed,rcpack);
        bin = hashkmer(packed,0);
        memcpy(kmerBuf[bin] + nwords*binTally[bin], packed, kmerSize);
        binTally[bin]++;
        if (binTally[bin] == maxKmersPerBin) serialize();
        kmer = rcpack;
    }
    bin = hashkmer(kmer,0);
    memcpy(kmerBuf[bin] + nwords*binTally[bin], kmer, kmerSize);
    binTally[bin]++;
    if (binTally[bin] == maxKmersPerBin) serialize();
    // pack the rest of the sequence
    for (size_t i=k; i<(size_t)length;i++) {
        nextKmer(packed, seq[i]);
        if (mode == CANONICAL)
            kmer = canonicalize(packed,rcpack);
        if (mode == BOTH) {
            kmer = canonicalize(packed,rcpack);
            bin = hashkmer(packed,0);
            memcpy(kmerBuf[bin] + nwords*binTally[bin], packed, kmerSize);
            binTally[bin]++;
            if (binTally[bin] == maxKmersPerBin) serialize();
            kmer = rcpack;
        }
        bin = hashkmer(kmer,0);
        memcpy(kmerBuf[bin] + nwords*binTally[bin], kmer, kmerSize);
        binTally[bin]++;
        if (binTally[bin] == maxKmersPerBin) serialize();
    }
}

void Kmerizer::save() {
    serialize();
    
    if (batches>1)
        mergeBatches();
    else
        for (size_t bin=0;bin<NBINS;bin++) {
            char ofname[100];
            char nfname[100];
            sprintf(ofname,"%s/%zi-mers.%zi.1",outdir,k,bin);
            sprintf(nfname,"%s/%zi-mers.%zi",outdir,k,bin);
            if (rename(ofname,nfname) != 0)
                perror("error renaming file");
            sprintf(ofname,"%s/%zi-mers.%zi.1.idx",outdir,k,bin);
            sprintf(nfname,"%s/%zi-mers.%zi.idx",outdir,k,bin);
            if (rename(ofname,nfname) != 0)
                perror("error renaming file");
        }
}

void Kmerizer::load() {
    fprintf(stderr,"Kmerizer::load()");
    clock_t start = clock();
    boost::thread_group tg;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doLoadIndex, this, i, j));
    }
    tg.join_all();
    state = QUERY;
    fprintf(stderr," took %f seconds\n",((float)(clock()-start))/CLOCKS_PER_SEC);
}

// given one kmer, pack it, canonicalize it, hash it, find it
uint32_t Kmerizer::find(const char* seq) {
    // pack it
    kword_t packed[nwords];
    kword_t rcpack[nwords];
    memset(packed,0,kmerSize);
    for (size_t i = 0; i < k; i++)
        nextKmer(packed,seq[i]);
    kword_t *kmer = packed;

    // canonicalize it
    if (mode == CANONICAL) {
        kmer = canonicalize(packed,rcpack);
    }

    // hash it
    size_t bin = hashkmer(kmer,0);

    // check if we've loaded the bit slices for this bin
    if (slices[bin].empty()) {
        // load the kmer slices for this bin
        char fname[100];
        sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);
        vector<uint32_t> junk;
        readBitmap(fname,junk,slices[bin]);
    }

    // search the index by iterative boolean operations
    BitVector *res = new BitVector(true); // first create an empty bitvector
    res->appendFill(true,slices[bin][0]->getSize()); // then make it all 1's
    for (size_t w=0;w<nwords;w++) {
        for (size_t b=0;b<64;b++) {
            if (kmer[w] & (1ULL << 63-b))
                *res &= *(slices[bin][w*nwords + b]);
            else
                *res &= *(slices[bin][w*nwords + b]->copyflip());
            if (res->cnt() == 0) return 0;
        }
    }
    // if we've reached this point, there should be one set bit in res
    if (res->cnt() != 1) {
        fprintf(stderr,"expected only 1, but cnt() returned %u\n",res->cnt());
        exit(1);
    }
    // need to figure out which bit
    // and lookup the associated count
    res->decompress();
    vector<uint32_t> hits = res->getWords();
    return frequency(bin,hits[0]);
}

uint32_t Kmerizer::frequency(size_t bin, uint32_t pos) {
    size_t key=1;
    while (key<kmerFreq[bin].size())
        if (counts[bin][key]->find(pos))
            key++;
        else
            break;
    return kmerFreq[bin][key-1];
}


void Kmerizer::histogram() {
    if (batches>1)
        save();

    if (state == READING)
        uniqify();

    // merge the NBINS counts bvecs
    uint32_t todo = NBINS;
    size_t offset[NBINS];
    size_t done[NBINS];
    memset(offset,0,sizeof(size_t)*NBINS);
    memset(done,0,sizeof(size_t)*NBINS);
    uint32_t key=1; // min kmer frequency
    while (todo > 0) {
        uint32_t val = 0;
        for (size_t i=0;i<NBINS;i++) {
            if (done[i]==0) {
                if (kmerFreq[i][offset[i]] == key) {
                    val += counts[i][offset[i]]->cnt();
                    if (offset[i]<kmerFreq[i].size()-1) // undo the range encoding <=
                        val -= counts[i][offset[i]+1]->cnt();
                    offset[i]++;
                    if (offset[i] == kmerFreq[i].size()) {
                        done[i]=1;
                        todo--;
                    }
                }
            }
        }
        if (val > 0) {
            printf("%u %u\n",key,val);
        }
        key++;
    }
}

void Kmerizer::serialize() {
    batches++;
    // unique the batch
    uniqify();
    // write to disk
    writeBatch();

    // zero the data
    memset(binTally,0,sizeof(uint32_t)*NBINS);
    for (size_t i=0;i<NBINS;i++) {
        for (size_t j=0; j<counts[i].size(); j++)
            delete counts[i][j]; // BitVector destructor
    }
    state = READING;
}

void Kmerizer::uniqify() {
    fprintf(stderr,"Kmerizer::uniqify()");
    clock_t start = clock();
    boost::thread_group tg;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doUnique, this, i, j));
    }
    tg.join_all();
    state = QUERY;
    fprintf(stderr," took %f seconds\n",((float)(clock()-start))/CLOCKS_PER_SEC);
}

void Kmerizer::writeBatch() {
    fprintf(stderr,"Kmerizer::writeBatch()");
    clock_t start = clock();
    boost::thread_group tg;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doUnique, this, i, j));
    }
    tg.join_all();
    fprintf(stderr," took %f seconds\n",((float)(clock()-start))/CLOCKS_PER_SEC);
}

void Kmerizer::mergeBatches() {
    fprintf(stderr,"Kmerizer::mergeBatches()");
    clock_t start=clock();
    boost::thread_group tg;
    for (size_t i=0;i<NBINS;i+=threadBins) {
        size_t j=(i+threadBins>NBINS) ? NBINS : i+threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doMergeBatches, this, i, j));
    }
    tg.join_all();
    batches=1;
    fprintf(stderr," took %f seconds\n",((float)(clock()-start))/CLOCKS_PER_SEC);
}


int compare_kmers1(const void *k1, const void *k2) {
    if (*(kword_t*)k1 > *(kword_t*)k2) return 1;
    if (*(kword_t*)k1 < *(kword_t*)k2) return -1;
    return 0;
}
int compare_kmers2(const void *k1, const void *k2) {
    for (size_t i=0;i<2;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers3(const void *k1, const void *k2) {
    for (size_t i=0;i<3;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers4(const void *k1, const void *k2) {
    for (size_t i=0;i<4;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers5(const void *k1, const void *k2) {
    for (size_t i=0;i<5;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers6(const void *k1, const void *k2) {
    for (size_t i=0;i<6;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers7(const void *k1, const void *k2) {
    for (size_t i=0;i<7;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers8(const void *k1, const void *k2) {
    for (size_t i=0;i<8;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}

int kmercmp(const void *k1, const void *k2, size_t nwords) {
    for (size_t i=0;i<nwords;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
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


void Kmerizer::doUnique(const size_t from, const size_t to) {
    for (size_t bin=from; bin<to; bin++) {
        if (1) {
        // stl sort
            if (nwords==1)
                sort((kmer1_t*)kmerBuf[bin], (kmer1_t*)(kmerBuf[bin]) + binTally[bin]);
            if (nwords==2)
                sort((kmer2_t*)kmerBuf[bin], (kmer2_t*)(kmerBuf[bin] + binTally[bin]));
            if (nwords==3)
                sort((kmer3_t*)kmerBuf[bin], (kmer3_t*)(kmerBuf[bin] + binTally[bin]));
            if (nwords==4)
                sort((kmer4_t*)kmerBuf[bin], (kmer4_t*)(kmerBuf[bin] + binTally[bin]));
            if (nwords==5)
                sort((kmer5_t*)kmerBuf[bin], (kmer5_t*)(kmerBuf[bin] + binTally[bin]));
            if (nwords==6)
                sort((kmer6_t*)kmerBuf[bin], (kmer6_t*)(kmerBuf[bin] + binTally[bin]));
            if (nwords==7)
                sort((kmer7_t*)kmerBuf[bin], (kmer7_t*)(kmerBuf[bin] + binTally[bin]));
            if (nwords==8)
                sort((kmer8_t*)kmerBuf[bin], (kmer8_t*)(kmerBuf[bin] + binTally[bin]));
        }
        if (0) {
            // qsort
            switch(nwords) {
                case 1:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers1);
                case 2:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers2);
                case 3:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers3);
                case 4:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers4);
                case 5:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers5);
                case 6:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers6);
                case 7:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers7);
                case 8:
                    qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers8);
            }
        }
        // uniq
        uint32_t distinct = 0;
        vector<uint32_t> tally;
        tally.push_back(1); // first kmer
        for (size_t i=1;i<binTally[bin];i++) {
            kword_t *ith = kmerBuf[bin] + i*nwords;
            if(kmercmp(kmerBuf[bin] + distinct*nwords, ith, nwords) == 0)
                tally.back()++;
            else {
                distinct++;
                tally.push_back(1);
                memcpy(kmerBuf[bin] + distinct*nwords, ith, kmerSize);
            }
        }
//      fprintf(stderr,"uniqify[%zi] reduced from %u - %u = %u\n",bin,binTally[bin],distinct+1,binTally[bin] - distinct - 1);
        binTally[bin] = distinct+1;
        // create a bitmap index for the tally vector
        rangeIndex(tally, kmerFreq[bin], counts[bin]);
    }
}

void Kmerizer::printKmer(kword_t * kmer) {
    size_t bpw = 8 * sizeof(kword_t);
    for (size_t j = 0; j < nwords; j++) {
        fprintf(stderr," ");
        for (size_t b=0;b<bpw;b++) {
            fprintf(stderr,"%d",((*(kmer+j) & (1ULL << (bpw-b-1)))) ? 1 : 0);
        }
    }
    fprintf(stderr,"\n");
}

void
Kmerizer::bitSlice(
    kword_t *    kmers,
    const size_t n,
    BitVector **      kmer_slices,
    size_t       nbits)
{
    // initialize WAH compressed bitvectors
    for (size_t i = 0; i < nbits; i++) {
        kmer_slices[i] = new BitVector(true);
    }
    bitset<256> bbit; // for keeping track of the set bits
    size_t boff[nbits]; // beginning of the current run of 1's or 0's
    memset(boff,0,nbits*sizeof(size_t)); // initialize to 0
    unsigned int bpw = 8*sizeof(kword_t); // bits per word

    // mark the set bits in the first kmer
    if (nwords>1) {
        for (size_t w=0;w<nwords;w++) {
            unsigned int count = popCount(kmers[w]);
            for (unsigned int r=1; r<=count; r++)
                bbit.set(selectBit(kmers[w],r) + w*bpw,1);
        }
        for (size_t i=1;i<n;i++) {
            // when n is large the number of different bits between kmer i and kmer i-1 is small
            // so use xor, popCount, and selectBit to identify the changed bit positions
            kword_t * kmer = kmers + i * nwords;
            kword_t * prev = kmer - nwords;
            for (size_t w=0;w<nwords;w++) {
                kword_t x = kmer[w] ^ prev[w];
                unsigned int count = popCount(x);
                for (unsigned int r = 1; r<=count; r++) {
                    unsigned int b = selectBit(x,r) + w*bpw;
                    kmer_slices[b]->appendFill(bbit.test(b),i-boff[b]);
                    bbit.flip(b);
                    boff[b] = i;
                }
            }
        }
    }
    else { // stripped down version for nwords=1
        unsigned int count = popCount(kmers[0]);
        for (unsigned int r=1; r<=count; r++)
            bbit.set(selectBit(kmers[0],r),1);
        for (size_t i=1;i<n;i++) {
            // when n is large the number of different bits between kmer i and kmer i-1 is small
            // so use xor, popCount, and selectBit to identify the changed bit positions
            uint32_t x = kmers[i] ^ kmers[i-1];
            unsigned int count = popCount(x);
            for (unsigned int r = 1; r<=count; r++) {
                unsigned int b = selectBit(x,r);
                kmer_slices[b]->appendFill(bbit.test(b),i-boff[b]);
                bbit.flip(b);
                boff[b] = i;
            }
        }
    }
    // append the last runs
    for (size_t b=0; b<nbits; b++)
        if(boff[b] < n-1)
            kmer_slices[b]->appendFill(bbit.test(b),n-boff[b]);
}

void Kmerizer::doWriteBatch(const size_t from, const size_t to) {
    for (size_t bin=from; bin<to; bin++) {
        // convert kmerBuf[bin] into a bit-sliced bitmap index
//      fprintf(stderr,"binTally[%zi]: %u\n",bin,binTally[bin]);
        const size_t nbits = 8*kmerSize;
        BitVector* kmer_slices[nbits];
        bitSlice(kmerBuf[bin],binTally[bin],kmer_slices,nbits);
        
        FILE *fp;
        if (0) {
            // open output file for kmerBuf
            char kmer_file[100];
            sprintf(kmer_file,"%s/%zi-mers.%zi.%zi",outdir,k,bin,batches);
            fp = fopen(kmer_file, "wb");
            fwrite(kmerBuf[bin],kmerSize,binTally[bin],fp);
            fclose(fp);
        }
        if (1) {
            // write the kmer_slices to a file
            char kmer_file[100];
            sprintf(kmer_file,"%s/%zi-mers.%zi.%zi",outdir,k,bin,batches);
            fp = fopen(kmer_file, "wb");
            size_t n_slices = 8*sizeof(kword_t);
            fwrite(&n_slices,sizeof(size_t),1,fp);
            for (size_t b=0;b<n_slices;b++) {
                uint32_t c = kmer_slices[b]->cnt();
                fwrite(&c,sizeof(uint32_t),1,fp);
            }
            // fprintf(stderr,"binTally[%zi]: %u\n",bin,binTally[bin]);
            for (size_t b=0;b<n_slices;b++) {
                uint32_t *buf;
                size_t bytes = kmer_slices[b]->dump(&buf);
                // fprintf(stderr,"%zi: %zi,",b,bytes);
                fwrite(&bytes,sizeof(size_t),1,fp);
                fwrite(buf,1,bytes,fp);
                free(buf);
            }
//          fprintf(stderr,"\n");
            fclose(fp);
        }
        // open output file for counts
        char counts_file[100];
        sprintf(counts_file,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,batches);
        fp = fopen(counts_file, "wb");
        // first write the number of distinct values
        // then write the distinct values
        // for each distinct value, write out the BitVector (size,count,rle,words.size(),words)
        size_t n_distinct = kmerFreq[bin].size();
        fwrite(&n_distinct,sizeof(size_t),1,fp);
        fwrite(kmerFreq[bin].data(),sizeof(uint32_t),n_distinct,fp);
        for (size_t i=0;i<n_distinct;i++) {
            uint32_t *buf;
            size_t bytes = counts[bin][i]->dump(&buf);
            // how many bytes are we writing
            fwrite(&bytes,sizeof(size_t),1,fp);
            // write them
            fwrite(buf,1,bytes,fp);
            free(buf);
        }
        fclose(fp);
    }
}

// when reading counts, values is the array of kmerFrequencies
// when reading kmers, values is the number of set bits in each bit slice
void Kmerizer::readBitmap(const char* idxfile, vector<uint32_t> &values, vector<BitVector*> &index) {
    // open idxfile
    FILE *fp;
    fp = fopen(idxfile,"rb");
    // fread to repopulate values and bvecs
    size_t n_distinct;
    fread(&n_distinct,sizeof(size_t),1,fp);
    values.resize(n_distinct);
    fread(values.data(),sizeof(uint32_t),n_distinct,fp); // assumes 32-bit words in bvec
    index.resize(n_distinct);
    for (size_t i=0;i<n_distinct;i++) {
        size_t bytes;
        fread(&bytes,sizeof(size_t),1,fp);
        size_t words = bytes/sizeof(uint32_t);
        uint32_t buf[words];
        fread(buf,1,bytes,fp);
        index[i] = new BitVector(buf);
    }
    fclose(fp);
}

void Kmerizer::doLoadIndex(const size_t from, const size_t to) {
    for (size_t bin=from; bin<to;bin++) {
        char fname[100];
        sprintf(fname,"%s/%zi-mers.%zi.idx",outdir,k,bin);
        readBitmap(fname,kmerFreq[bin],counts[bin]);
    }
}

void Kmerizer::doMergeBatches(const size_t from, const size_t to) {
    for (size_t bin=from; bin<to; bin++) {
        // read the counts and kmers for each batch
        vector<BitVector*>    batch_counts[batches];
        vector<uint32_t> batch_values[batches];
        vector<BitVector*>    batch_slices[batches];
        vector<uint32_t> batch_slice_cnts[batches];
        kword_t  kmers[batches*nwords]; // next kmer in each active batch
        uint32_t btally[batches]; // frequency of next kmer in each batch
        size_t   offset[batches]; // keep track of position in each batch

        for (size_t i=0;i<batches;i++) {
            char fname[100];
            sprintf(fname,"%s/%zi-mers.%zi.%zi",outdir,k,bin,i+1);
            readBitmap(fname,batch_slice_cnts[i],batch_slices[i]);
            sprintf(fname,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,i+1);
            readBitmap(fname,batch_values[i],batch_counts[i]);
        }

        memset(offset, 0, sizeof(size_t)*batches);
        uint32_t todo = batches; // number of batches to process
        vector<uint32_t> tally;

        const size_t nbits = 8 * kmerSize;
        BitVector* merged_slices[nbits];
        bitset<256> bbit;
        unsigned int boff[nbits];
        unsigned int n=0;
        memset(boff, 0, nbits * sizeof(size_t));
        const unsigned int bpw = 8 * sizeof(kword_t); // bits per word

        for (size_t b=0;b<nbits;b++)
            merged_slices[b] = new BitVector(true);

        // read the first kmer and counts from each batch
        for (size_t i=0;i<batches;i++) {
            size_t err =
                pos2kmer(offset[i], kmers + i * nwords, batch_slices[i]);
            if (err != 0) {
                todo--;
                btally[i]=0;
            }
            else {
                btally[i] = pos2value(offset[i]++,batch_values[i],batch_counts[i]);
            }
        }
        // choose min
        size_t mindex = findMin(kmers,btally);
        uint32_t distinct[nwords];
        memcpy(distinct,kmers + mindex*nwords, kmerSize);
        tally.push_back(btally[mindex]);
        // mark the set bits in the first distinct kmer
        for (size_t w=0;w<nwords;w++) {
            unsigned int count = popCount(distinct[w]);
            for (unsigned int r=1; r<=count; r++)
                bbit.set(selectBit(distinct[w],r)+w*bpw,1);
        }
        // replace min
        // iterate until there's nothing left to do
        
        while(todo>0) {
            size_t err = pos2kmer(offset[mindex],kmers + mindex*nwords, batch_slices[mindex]);
            if (err != 0) {
                todo--;
                btally[mindex]=0;
                if (todo==0) break;
            }
            else
                btally[mindex] = pos2value(offset[mindex]++,batch_values[mindex],batch_counts[mindex]);
            mindex = findMin(kmers,btally);
            // compare to distinct
            if (kmercmp(kmers + mindex*nwords, distinct, nwords) == 0) // same kmer
                tally.back() += btally[mindex];
            else { // find the changed bits
                tally.push_back(btally[mindex]);
                n++;
                kword_t *minkmer = kmers + mindex*nwords;
                for (size_t w=0;w<nwords;w++) {
                    kword_t x = distinct[w] ^ minkmer[w];
                    unsigned int count = popCount(x);
                    for (unsigned int r = 1; r<=count; r++) {
                        unsigned int b = selectBit(x,r) + w*bpw;
                        merged_slices[b]->appendFill(bbit.test(b),n-boff[b]);
                        bbit.flip(b);
                        boff[b]=n;
                    }
                    distinct[w] = minkmer[w];
                }
            }
        }
        // finish the bitvectors
        n++;
        for (size_t b=0; b<nbits; b++)
            if(boff[b]<n-1)
                merged_slices[b]->appendFill(bbit.test(b),n-boff[b]);
        
        rangeIndex(tally,kmerFreq[bin],counts[bin]);

        // open an output file for the merged distinct kmers
        FILE *ofp;
        char fname[100];
        sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);
        ofp = fopen(fname,"wb");

        fwrite(&nbits,sizeof(size_t),1,ofp);
        for (size_t b=0;b<nbits;b++) {
            uint32_t c = merged_slices[b]->cnt();
            fwrite(&c,sizeof(uint32_t),1,ofp);
        }
        for (size_t b=0;b<nbits;b++) {
            uint32_t *buf;
            size_t bytes = merged_slices[b]->dump(&buf);
            // fprintf(stderr,"%zi: %zi,",b,bytes);
            fwrite(&bytes,sizeof(size_t),1,ofp);
            fwrite(buf,1,bytes,ofp);
            free(buf);
        }

        fclose(ofp);
        
        char counts_file[100];
        sprintf(counts_file,"%s/%zi-mers.%zi.idx",outdir,k,bin);
        ofp = fopen(counts_file, "wb");
        // first write the number of distinct values
        // then write the distinct values
        // for each distinct value, write out the BitVector (size,count,rle,words.size(),words)
        int n_distinct = kmerFreq[bin].size();
        fwrite(&n_distinct,sizeof(int),1,ofp);
        fwrite(kmerFreq[bin].data(),sizeof(uint32_t),n_distinct,ofp);
        for (size_t i=0;i<n_distinct;i++) {
            uint32_t *buf;
            size_t bytes = counts[bin][i]->dump(&buf);
            fwrite(&bytes,sizeof(size_t),1,ofp);
            fwrite(counts[bin][i]->getWords().data(),1,counts[bin][i]->bytes(),ofp);
        }
        fclose(ofp);
        // delete input files
        // for (size_t i=0;i<batches;i++) {
        //  sprintf(fname,"%s/%zi-mers.%zi.%zi",outdir,k,bin,i+1);
        //  if (remove(fname) != 0) perror("error deleting file");
        //  sprintf(fname,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,i+1);
        //  if (remove(fname) != 0) perror("error deleting file");
        // }
    }
}

inline size_t
Kmerizer::findMin(const kword_t * kmers, const uint32_t * kcounts) {
    size_t mindex = 0;
    while (kcounts[mindex] == 0) mindex++;
    for (size_t i=mindex+1;i<batches;i++)
        if (kcounts[i] != 0)
            if (kmercmp(kmers + i*nwords,kmers + mindex*nwords,nwords) < 0)
                mindex = i;

    return mindex;
}

// for each distinct value in the vec create a bitvector
// indexing the positions in the vec holding a value <= v
void Kmerizer::rangeIndex(vector<uint32_t> &vec, vector<uint32_t> &values, vector<BitVector*> &index) {
    // find the distinct values in vec
    values.clear();
    // use a bitset to mark the distinct values (from 0-255)
    // and another vector for overflow (which we'll sort/uniq later)
    bitset<256> mybits;
    vector<uint32_t> overflow;
    vector<uint32_t>::iterator it;
    for (it = vec.begin(); it < vec.end(); ++it) {
        if (*it >= 256) 
            overflow.push_back(*it);
        else
            mybits[*it]=1;
    }
    // populate values vector with set bits in mybits
    for (uint32_t i=0;i<256;i++)
        if (mybits.test(i))
            values.push_back(i);
    
    if (overflow.size() > 0) {
        sort(overflow.begin(),overflow.end());
        it = unique(overflow.begin(),overflow.end());
        values.insert(values.end(),overflow.begin(),it);
    }

    // setup a vector for each range
    vector<uint32_t> vrange[values.size()];
    // iterate over the vec and push the offset onto each range
    for (size_t i=0;i<vec.size();i++) {
        it = lower_bound(values.begin(),values.end(),vec[i]);
        for (size_t j=0;j<=(it - values.begin());j++)
            vrange[j].push_back(i);
    }
    // create a bitvector for each range
    index.resize(values.size());
    for (size_t i=0;i<values.size();i++)
        index[i] = new BitVector(vrange[i]);
}

// need an iterator that works with a BitVector mask BitVector->next_one(); - returns position of set bit
// or number of bits (invalid response)

// We can do this in parallel if we figure out in advance the offset within the output file for each bin.
// A frequency histogram can get us the number of 1 digit counts, 2 digit counts, 3 digit counts, etc.
// So, if there are mask->cnt() set bits, we need n*(k+2) + n1 + 2n2 + 3n3 etc bytes. Each kmer takes k bytes, and you have a tab and a newline character on each line. n1 = 1 digit counts, n2 = 2 digit counts, etc...
void Kmerizer::dump(char *fname) {
    // create a bitvector of all 1's for each bin
    // and call dump(FILE *fp, BitVector *mask)
    BitVector *mask[NBINS];
    for (size_t i=0;i<NBINS;i++) {
        mask[i] = new BitVector(true); // first create an empty bitvector
        mask[i]->appendFill(true,counts[i][0]->getSize()); // then make it all 1's
    }
    sdump(fname,mask);
}

// find kmers with frequencies in the given range[min,max]
// when max < min, ignore max and find kmers with frequencies in the range[min,infinity]
void Kmerizer::filter(uint32_t min, uint32_t max, BitVector **mask) {
    // do range queries over each bin to build up an array of BitVector masks.
    boost::thread_group tg;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doFilter, this, i, j, min, max, mask));
    }
    tg.join_all();
}


void Kmerizer::pdump(char *fname, BitVector **mask) {
    // determine the uncompressed bin sizes in advance so we can doDump in parallel
    // binsize[bin] = mask[bin]->cnt() * (k+2) + d*number of d digit counts
    uint32_t binsize[NBINS];
    uint32_t total_size=0;
    for (size_t bin=0;bin<NBINS;bin++) {
        binsize[bin] = mask[bin]->cnt() * (k+2); // k-mer\t\n
        // d digit counts
        uint32_t d=1;
        vector<uint32_t>::iterator it = lower_bound(kmerFreq[bin].begin(),kmerFreq[bin].end(),10);
        while(it != kmerFreq[bin].end()) {
            if (it > kmerFreq[bin].begin()) {
                it--;
                BitVector *dmask = *counts[bin][*it] & *mask[bin];
                binsize[bin] += d*dmask->cnt();
            }
            d++;
            it = lower_bound(it,kmerFreq[bin].end(),10^d);
        }
        total_size += binsize[bin];
    }
    // use memory mapped IO to do this
    // maybe open the output file, resize it, and then mmap it to some buffer
    // You can't do this with stdout. is that a problem?
    char *buff;
//  mmap(buff, total_size);
    boost::thread_group tg;
    uint32_t offset=0;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        // pass a pointer to the right place in the mmap'd output file
        tg.create_thread(boost::bind(&Kmerizer::doPdump, this, i, j, buff + offset, mask));
        for (size_t bin=i;bin<j;bin++)
            offset += binsize[bin];
    }
    tg.join_all();
    // close the output file
}

void Kmerizer::sdump(char *fname, BitVector **mask) {
    // open output file
    FILE *fp;
    fp = fopen(fname, "w");
    size_t bpw = 8*sizeof(kword_t); // bits per word
    char kstr[k+1]; // unpack each kmer into this char array.
    for (size_t bin=0;bin<NBINS;bin++) {
        // check if bit slices have been loaded
        if (slices[bin].empty()) {
            char bfname[100];
            sprintf(bfname,"%s/%zi-mers.%zi",outdir,k,bin);
            vector<uint32_t> tmp;
            readBitmap(bfname,tmp,slices[bin]);
        }
        // iterate over the set bits in mask[bin]
        uint32_t bv_len = mask[bin]->getSize();
        uint32_t next1 = mask[bin]->nextOne(0);
        while (next1 < bv_len) {
            kword_t kmer[nwords];
            for (size_t w=0;w<nwords;w++) {
                kmer[w]=0;
                for (size_t b=0;b<bpw;b++) {
                    if(slices[bin][w*bpw + b]->find(next1)) {
                        kmer[w] |= 1 << bpw-b-1;
                    }
                }
            }
            // output the kmer and count
            unpack(kmer, kstr);
            fprintf(fp,"%s\t%u\n",kstr,frequency(bin,next1));
            uint32_t next1 = mask[bin]->nextOne(next1+1);
        }
    }
    fclose(fp);
}

void
Kmerizer::doPdump(const size_t from,
                 const size_t to,
                 char *       buff,
                 BitVector ** mask)
{
    size_t bpw = 8 * sizeof(kword_t); // bits per word
    char kstr[k+1]; // unpack each kmer into this char array.
    for (size_t bin = from; bin < to; bin++) {
        // check if we've loaded the bit slices for this bin
        if (slices[bin].empty()) {
            // load the kmer slices for this bin
            char fname[100];
            sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);
            vector<uint32_t> junk;
            readBitmap(fname,junk,slices[bin]);
        }
        // iterate over the set bits in mask[bin]
        uint32_t bv_len = mask[bin]->getSize();
        uint32_t next1 = mask[bin]->nextOne(0);
        while (next1 < bv_len) {
            kword_t kmer[nwords];
            for (size_t w=0;w<nwords;w++) {
                kmer[w]=0;
                for (size_t b=0;b<bpw;b++) {
                    if(slices[bin][w*bpw + b]->find(next1)) {
                        kmer[w] |= 1 << bpw-b-1;
                    }
                }
            }
            // output the kmer and count
            unpack(kmer, kstr);
            // instead of fprintf, use memcpy() to write kstr and frequency to the mmap'd output file
//          fprintf(fp,"%s\t%u\n",kstr,frequency(bin,next1));
            uint32_t next1 = mask[bin]->nextOne(next1+1);
        }
    }
}

// range encoding: bitvector x marks kmers that occur <= x times
// therefore counts[bin].back() marks all rows
void
Kmerizer::doFilter(const size_t from,
                   const size_t to,
                   uint32_t min,
                   uint32_t max,
                   BitVector **mask)
{
    for (size_t bin = from; bin < to; bin++) {
        if (max > 0) {
            vector<uint32_t>::iterator it =
                lower_bound(kmerFreq[bin].begin(),kmerFreq[bin].end(),max);
            if (it == kmerFreq[bin].begin()) {
                // create a null mask and set the length
                mask[bin] = new BitVector(true);
                mask[bin]->appendFill(false,counts[bin][0]->getSize());
            }
            else {
                if (*it < max) it--;
                mask[bin]->copy(counts[bin][it - kmerFreq[bin].begin()]);
            }
        }
        else
            mask[bin]->copy(counts[bin].back()); // this is how we deal with max <= 0 (no max)
        if (mask[bin]->cnt() > 0 && min>1) { // nand
            // need to find the first BitVector thats < min
            vector<uint32_t>::iterator it = lower_bound(kmerFreq[bin].begin(),kmerFreq[bin].end(),min);
            if (it > kmerFreq[bin].begin()) {
                it--;
                *mask[bin] &= *(counts[bin][it - kmerFreq[bin].begin()]->copyflip());
            }
        }
    }
}
