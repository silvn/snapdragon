#ifndef SNAPDRAGON_KMERIZER_H
#define SNAPDRAGON_KMERIZER_H

#define DEBUG false
#define NBINS 256
#define CANONICAL 'C'
#define BOTH 'B'
#define COUNTING 'C'
#define SAVING 'S'

#include <vector>
#include "../boost/threadpool.hpp"

typedef uint64_t kword_t;
using namespace std;

class Kmerizer {
    size_t  k;
    kword_t kmask;
    size_t  shiftlastby;
    size_t  rshift, lshift;
    size_t  nwords;
    size_t  kmerSize; // in bytes
    size_t  nthreads;
    size_t  maxSeqLength; // for parallelization
    char    mode;
    char    state;
    char *  outdir;

    // lookup tables for common kmers
    kword_t  * kmerLutK[NBINS]; // key = kmer
    uint32_t * kmerLutV[NBINS]; // value = frequency
    uint32_t   kmerLutS[NBINS]; // size of lookup table

    // buffers for rare kmers
    kword_t  * kmerBuf[NBINS];
    uint32_t   kmerBufTally[NBINS]; // number of kmers stored
    uint32_t   kmerBufSize[NBINS];  // capacity of buffer

    uint32_t   batches[NBINS]; // keep track of the number of serialized batches
    
    // some semaphores to protect the kmer lookup tables and buffers
    int usingLut[NBINS]; // number of threads using the lookup table (-1 means don't edit)
    int usingBuf[NBINS]; // number of threads using the buffer

    boost::threadpool::pool tp;

public:
    // constructor - sets up worker threads
    Kmerizer(const size_t klength,
             const size_t nthreads,
             const char * outdir,
             const char   mode);

    ~Kmerizer() {};

    // allocate memory for counting
    int allocate(const size_t maximem);

    // extract kmers from the sequence
    // seq should only contain non-ACGT characters
    void addSequence(const char* seq, const int length);

    // write distinct kmers and counts to disk (merging multiple batches)
    void save();

    void load() {}
    void histogram() {}

private:

    // worker thread function to extract (canonical) kmers from a sequence of [ACGT]
    // this function relies on addSequence() to split sequences into managable chunks
    // (length <= maxSeqLength)
    void addSeq(const char* seq, const int length);
    // version for when nwords==1
    void addSeq1(const char* seq, const int length);

    void insertKmer(size_t bin, kword_t *kmer);
    void insertKmer1(size_t bin, kword_t *kmer);
    
    int searchLut(kword_t *kmer, size_t bin);
    int searchLut1(kword_t *kmer, size_t bin);

    // update kmer given the next nucl by shifting by two bits
    size_t nextKmer(kword_t* kmer, size_t bin, const char nucl);
    size_t nextKmer1(kword_t* kmer, size_t bin, const char nucl);

    // choose the lesser between a kmer and its reverse complement
    kword_t* canonicalize(kword_t *packed, kword_t *rcpack, size_t *bin) const;
    kword_t* canonicalize1(kword_t *packed, kword_t *rcpack, size_t *bin) const;

    // pack nucleotides into 2 bits
    inline kword_t twoBit(const unsigned char nuc) const;
    
    // reverse complement one kword_t
    inline kword_t revcomp(const kword_t val) const;
    // only do the reverse complement of the 8 least significant digits
    inline kword_t revcomp_8(const kword_t val) const;
    
    void uniqify(size_t bin);
    void updateLut(size_t bin, vector<uint32_t> &tally);
    
    void saveBin(size_t bin);
};

inline kword_t Kmerizer::twoBit(const unsigned char nuc) const {
    static const kword_t table[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    return table[nuc];
}

inline kword_t Kmerizer::revcomp(const kword_t val) const {
    // reverse complement lookup table (8 bits at a time)
    static const kword_t rctable[256] =
    {
        255,191,127,63,239,175,111,47,223,159,95,31,207,143,79,15,
        251,187,123,59,235,171,107,43,219,155,91,27,203,139,75,11,
        247,183,119,55,231,167,103,39,215,151,87,23,199,135,71,7,
        243,179,115,51,227,163,99,35,211,147,83,19,195,131,67,3,
        254,190,126,62,238,174,110,46,222,158,94,30,206,142,78,14,
        250,186,122,58,234,170,106,42,218,154,90,26,202,138,74,10,
        246,182,118,54,230,166,102,38,214,150,86,22,198,134,70,6,
        242,178,114,50,226,162,98,34,210,146,82,18,194,130,66,2,
        253,189,125,61,237,173,109,45,221,157,93,29,205,141,77,13,
        249,185,121,57,233,169,105,41,217,153,89,25,201,137,73,9,
        245,181,117,53,229,165,101,37,213,149,85,21,197,133,69,5,
        241,177,113,49,225,161,97,33,209,145,81,17,193,129,65,1,
        252,188,124,60,236,172,108,44,220,156,92,28,204,140,76,12,
        248,184,120,56,232,168,104,40,216,152,88,24,200,136,72,8,
        244,180,116,52,228,164,100,36,212,148,84,20,196,132,68,4,
        240,176,112,48,224,160,96,32,208,144,80,16,192,128,64,0,
    };
    
    return
        (rctable[val&0xFFUL]<<56) |
        (rctable[(val>>8)&0xFFUL]<<48) |
        (rctable[(val>>16)&0xFFUL]<<40) |
        (rctable[(val>>24)&0xFFUL]<<32) |
        (rctable[(val>>32)&0xFFUL]<<24) |
        (rctable[(val>>40)&0xFFUL]<<16) |
        (rctable[(val>>48)&0xFFUL]<<8) |
        (rctable[(val>>56)&0xFFUL]);
}

// just return the reverse complement of the low 8 bits
inline kword_t Kmerizer::revcomp_8(const kword_t val) const {
    // reverse complement lookup table (8 bits at a time)
    static const kword_t rctable[256] =
    {
        255,191,127,63,239,175,111,47,223,159,95,31,207,143,79,15,
        251,187,123,59,235,171,107,43,219,155,91,27,203,139,75,11,
        247,183,119,55,231,167,103,39,215,151,87,23,199,135,71,7,
        243,179,115,51,227,163,99,35,211,147,83,19,195,131,67,3,
        254,190,126,62,238,174,110,46,222,158,94,30,206,142,78,14,
        250,186,122,58,234,170,106,42,218,154,90,26,202,138,74,10,
        246,182,118,54,230,166,102,38,214,150,86,22,198,134,70,6,
        242,178,114,50,226,162,98,34,210,146,82,18,194,130,66,2,
        253,189,125,61,237,173,109,45,221,157,93,29,205,141,77,13,
        249,185,121,57,233,169,105,41,217,153,89,25,201,137,73,9,
        245,181,117,53,229,165,101,37,213,149,85,21,197,133,69,5,
        241,177,113,49,225,161,97,33,209,145,81,17,193,129,65,1,
        252,188,124,60,236,172,108,44,220,156,92,28,204,140,76,12,
        248,184,120,56,232,168,104,40,216,152,88,24,200,136,72,8,
        244,180,116,52,228,164,100,36,212,148,84,20,196,132,68,4,
        240,176,112,48,224,160,96,32,208,144,80,16,192,128,64,0,
    };
    
    return rctable[val&0xFFUL];
}

#endif // #ifndef SNAPDRAGON_KMERIZER_H
