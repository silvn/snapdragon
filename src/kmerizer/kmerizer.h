#ifndef SNAPDRAGON_KMERIZER_H
#define SNAPDRAGON_KMERIZER_H

#define DEBUG false
#define NBINS 256
#define CANONICAL 'C'
#define BOTH 'B'
#define ASIS 'A'
#define COUNT 'C'
#define SAVE 'S'
#define MERGE 'M'
#define QUERY 'Q'
#include <queue>
#include <utility>
#include "../boost/threadpool.hpp"
#include "../bitmap/bitSlicedIndex.hpp"
#include "../bitmap/lcbitSlicedIndex.hpp"
// #include "../bitmap/rangeEncodedIndex.hpp"

typedef uint64_t kword_t;
using namespace std;
struct kmer1_t { kword_t a; };
struct kmer2_t { kword_t a[2]; };
struct kmer3_t { kword_t a[3]; };
struct kmer4_t { kword_t a[4]; };
struct kmer5_t { kword_t a[5]; };
struct kmer6_t { kword_t a[6]; };
struct kmer7_t { kword_t a[7]; };
struct kmer8_t { kword_t a[8]; };
struct kmer9_t { kword_t a[9]; };

class Kmerizer {
    size_t  k;
    kword_t kmask;
    size_t  shiftlastby;
    size_t  rshift, lshift;
    size_t  nwords;
    size_t  kmerSize; // in bytes
    size_t  nthreads;
    char    mode;
    char state; // see #define (COUNT, SAVE, QUERY, etc.)
    char *  outdir;
    bool    uniformBins;

    // buffers for kmers
    kword_t * kmerBuf [NBINS];
    vector<kword_t> kmerBuf1[NBINS];
    uint32_t   kmerBufTally[NBINS]; // number of kmers stored
    uint32_t   kmerBufSize[NBINS];  // capacity of buffer

    uint32_t   batches[NBINS]; // keep track of the number of serialized batches

    boost::threadpool::pool tp;
    
    BitVector<uint32_t> *BitMask[NBINS]; // query results (default all)
    char filtered[NBINS];
public:
    // constructor
    Kmerizer(const size_t klength,
             const size_t nthreads,
             const char * outdir,
             const char   mode);

    ~Kmerizer() {};

    // allocate memory for counting, setup ProducerConsumerQueues and schedule consumeKmers() jobs
    int allocate(const size_t maximem);

    // extract kmers from the sequence
    // seq should only contain non-ACGT characters
    void addSequence(const char* seq, const int length);

    // write distinct kmers and counts to disk (merging multiple batches)
    void save();
    
    // called after filter()
    void saveSubset(char *fname);

    // query functions
    void filter(char *colname, uint32_t min, uint32_t max); // results are stored in BitMask
    void query(const char* seq, const int length); // extracts kmers, queries indexes, stores results in BitMask
    void intersect(const char * lhs, const char * rhs); // for MANY kmers, intersect databases

    // BitMask is used by all of these functions
    void stats(char *colname); // unique, distinct, total, max_count
    void histo(char *colname); // for each frequency, report the number of kmers
    void dump(vector<char *> &columns); // write tab delimited text or fasta to stdout

private:

    // function to extract (canonical) kmers from a sequence of [ACGT]
    void addSeq(const char* seq, const int length);
    void addSeq1(const char* seq, const int length);

    void unpack(kword_t *kmer, size_t bin, char *kmerStr);

    // k-mers are not distributed evenly among the bins
    // redistributes free space among the bins
    void resizeBins1();
    uint32_t reservedSpace;

    // update kmer given the next nucl by shifting by two bits, returns nextKmer's bin
    size_t nextKmer (kword_t* kmer, size_t bin, const char nucl);
    size_t nextKmer1(kword_t* kmer, size_t bin, const char nucl);

    // choose the lesser between a kmer and its reverse complement
    kword_t* canonicalize (kword_t *packed, kword_t *rcpack, size_t *bin);
    kword_t* canonicalize1(kword_t *packed, kword_t *rcpack, size_t *bin);
    
    // search lookup table (kmerLutK[bin]) and increment count (kmerLutV[bin])
    // or save kmer for later (kmerBuf[bin])
    void insertKmer(size_t bin, kword_t *kmer);
    void insertKmer1(size_t bin, kword_t *kmer);
        
    // sort and count the kmers in kmerBuf[bin], updateLut (optional), write distinct kmers and counts to disk
    void uniqify(size_t bin);
    void uniqify1(size_t bin);

    // sort and count the kbt kmers in kb and put counts in tally
    uint32_t sortCount(kword_t *kb, uint32_t kbt, vector<uint32_t> &tally);
    // specialization for nwords==1
    uint32_t sortCount1(kword_t *kb, uint32_t kbt, vector<uint32_t> &tally);
    void sortCount1(vector<kword_t> &kb, vector<uint32_t> &tally);
    // partition kmers into bins before calling sortCount1()
    uint32_t binSortCount1(kword_t *kb, uint32_t kbt, vector<uint32_t> &tally);
    void binSortCount1(vector<kword_t> &kb, vector<uint32_t> &tally);
    void binSortCount1(vector<kword_t> &kb, BitSlicedIndex<kword_t> *kmerIdx, LCBitSlicedIndex<uint32_t> *countIdx);
    void sortCount1(vector<kword_t> &kb, BitSlicedIndex<kword_t> *kmerIdx, LCBitSlicedIndex<uint32_t> *countIdx);

    void writeBatch(size_t bin, vector<uint32_t> &tally);
    void sortwriteBatch(size_t bin);
    
    // call uniqify again and write lookup table to disk
    void saveBin(size_t bin);
    
    // merge raw kmers and counts into bitmap self-indexes
    void mergeBin(size_t bin);
    void mergeBin1(size_t bin);
    void filterBin(size_t bin, char *colname, uint32_t min, uint32_t max);
    static uint64_t reverse_complement(uint64_t v) {
      v = ((v >> 2)  & 0x3333333333333333UL) | ((v & 0x3333333333333333UL) << 2);
      v = ((v >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((v & 0x0F0F0F0F0F0F0F0FUL) << 4);
      v = ((v >> 8)  & 0x00FF00FF00FF00FFUL) | ((v & 0x00FF00FF00FF00FFUL) << 8);
      v = ((v >> 16) & 0x0000FFFF0000FFFFUL) | ((v & 0x0000FFFF0000FFFFUL) << 16);
      v = ( v >> 32                        ) | ( v                         << 32);
      return v ^ 0xFFFFFFFFFFFFFFFFUL;
    }

    // just return the reverse complement of the low 8 bits
    static uint64_t revcomp_8(uint64_t val) {
        // reverse complement lookup table (8 bits at a time)
        static const uint64_t rctable[256] =
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
    
        return rctable[val & 255UL];
    }
};
void binSortCounter(vector<kword_t> &kb, unsigned int my_k, vector<uint32_t> &tally);
void sortCountKmers(vector<kword_t> &kmers, vector<uint32_t> &tally);
#include "kmerizer.tpp"
#endif // #ifndef SNAPDRAGON_KMERIZER_H
