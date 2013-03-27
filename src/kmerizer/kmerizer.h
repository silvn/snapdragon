#ifndef KMERIZER_H
#define KMERIZER_H

#define DEBUG false
#define NBINS 256
#define CANONICAL 'C'
#define BOTH 'B'
#define READING 1
#define QUERY 2

#include <vector>
#include "bvec32.h"

typedef uint64_t word_t;
using namespace std;

class kmerizer {
	size_t k;
	word_t kmask;
	size_t shiftlastby;
	size_t rshift;
	size_t nwords;
	size_t kmer_size; // in bytes
	size_t threads;
	size_t thread_bins;
	size_t batches;
	size_t max_kmers_per_bin;
	char mode;
	char state;
	char *outdir;

	word_t* kmer_buf[NBINS]; // raw unsorted padded kmers, or sort|uniq'ed kmers
	uint32_t bin_tally[NBINS]; // number of kmers in each bin (or number of distinct kmers)
	vector<uint32_t> kmer_freq[NBINS]; // sorted distinct kmer frequencies
	vector<bvec32*> counts[NBINS]; // bitmap index of frequency counts

public:
	kmerizer(const size_t _k, const size_t _threads, char* _outdir, const char _mode);
	int allocate(const size_t maximem); // allocates memory for each kmer_buf
	void addSequence(const char* seq,const int length); // extract (canonicalized) kmers from the sequence
	void save(); // writes distinct kmers and rle counts to disk (merging multiple batches)
	void histogram(); // output the kmer count frequency distribution
	~kmerizer() {};

private:
	inline word_t twobit(const word_t val) const; // pack nucleotides into 2 bits
	inline void next_kmer(word_t* kmer, const char nucl); // shift kmer to make room for nucl
	inline void unpack(word_t* kmer, char *seq);
	inline word_t revcomp(const word_t val) const; // reverse complement
	inline uint8_t hashkmer(const word_t *kmer, const uint8_t seed) const; // to select a bin
	inline word_t* canonicalize(word_t *packed, word_t *rcpack) const;
	void serialize(); // kmer_buf is full. uniqify and write batch to disk
	void uniqify(); // qsort each kmer_buf, update bin_tally, and fill counts
	void do_unique(const size_t from, const size_t to); // for parallelization
	void writeBatch();
	void do_writeBatch(const size_t from, const size_t to); // for parallelization
	void mergeBatches();
	void do_mergeBatches(const size_t from, const size_t to);
	// is this too generic to go here?
	void range_index(vector<uint32_t> &vec, vector<uint32_t> &values, vector<bvec32*> &index);
	void read_index(const char* idxfile,vector<uint32_t> &values, vector<bvec32*> &index);
	size_t find_min(const word_t* kmers, const uint32_t* kcounts);
	uint32_t pos2value(size_t pos, vector<uint32_t> &values, vector<bvec32*> &index);
};

inline word_t kmerizer::twobit(const word_t val) const {
	static const uint8_t table[256] =
	{
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
	return (word_t)table[val];
}

inline void kmerizer::next_kmer(word_t* kmer, const char nucl) {
	// shift first kmer by 2 bits
	kmer[0] <<= 2;
	for(size_t w=1;w<nwords-1;w++) { // middle (full length) words
		kmer[w-1] |= kmer[w] >> 62; // move the top two bits
		kmer[w] <<= 2; // make room for the next word
	}
	// last word (if not the first)
	if(nwords>1) {
		kmer[nwords-2] |= kmer[nwords-1] >> shiftlastby;
		kmer[nwords-1] <<= 2;
	}
	kmer[nwords-1] |= twobit(nucl);
	kmer[nwords-1] &= kmask;
}

inline void kmerizer::unpack(word_t* kmer, char* seq) {
	static const char table[4] = {65, 67, 71, 84};
	for(size_t i=0;i<k;i++) {
		// which word, which nucl
		size_t w = i>>5;
		if (w == nwords-1) // not a full length word
			seq[i] = table[(kmer[w] >> (2*(k%32 - i%32) - 2)) & 3];
		else
			seq[i] = table[(kmer[w] >> (2*(i%32) - 2)) & 3];
	}
	seq[k] = '\0';
}

inline word_t kmerizer::revcomp(const word_t val) const {
	// bitwise reverse complement of values from 0 to 255
	static const uint8_t table[256] =
	{
		255,127,191,63,223,95,159,31,239,111,175,47,207,79,143,15,
		247,119,183,55,215,87,151,23,231,103,167,39,199,71,135,7,
		251,123,187,59,219,91,155,27,235,107,171,43,203,75,139,11,
		243,115,179,51,211,83,147,19,227,99,163,35,195,67,131,3,
		253,125,189,61,221,93,157,29,237,109,173,45,205,77,141,13,
		245,117,181,53,213,85,149,21,229,101,165,37,197,69,133,5,
		249,121,185,57,217,89,153,25,233,105,169,41,201,73,137,9,
		241,113,177,49,209,81,145,17,225,97,161,33,193,65,129,1,
		254,126,190,62,222,94,158,30,238,110,174,46,206,78,142,14,
		246,118,182,54,214,86,150,22,230,102,166,38,198,70,134,6,
		250,122,186,58,218,90,154,26,234,106,170,42,202,74,138,10,
		242,114,178,50,210,82,146,18,226,98,162,34,194,66,130,2,
		252,124,188,60,220,92,156,28,236,108,172,44,204,76,140,12,
		244,116,180,52,212,84,148,20,228,100,164,36,196,68,132,4,
		248,120,184,56,216,88,152,24,232,104,168,40,200,72,136,8,
		240,112,176,48,208,80,144,16,224,96,160,32,192,64,128,0
	};
	
	return
		((word_t)table[val&0xFFUL]<<56) |
		((word_t)table[(val>>8)&0xFFUL]<<48) |
		((word_t)table[(val>>16)&0xFFUL]<<40) |
		((word_t)table[(val>>24)&0xFFUL]<<32) |
		((word_t)table[(val>>32)&0xFFUL]<<24) |
		((word_t)table[(val>>40)&0xFFUL]<<16) |
		((word_t)table[(val>>48)&0xFFUL]<<8) |
		((word_t)table[(val>>56)&0xFFUL]);
}

inline uint8_t kmerizer::hashkmer(const word_t *kmer, const uint8_t seed) const {
	static const uint8_t Rand8[256] =
	{
		105,193,195, 26,208, 80, 38,156,128,  2,101,205, 75,116,139, 61,
		197,120,244, 51,185,132, 55,150,177,241,103,196, 13,237,136, 24,
		211, 56,207,  9, 30,145, 18,167,108, 32,106,151,  3, 54,248, 65,
		 17,198, 85, 95, 29, 83,  4,206,188,186,107,255,129, 35,142, 91,
		203,158, 74,138,162,135,102,114, 81,170, 10, 19,215, 57,214, 70,
		 37,163,231,227,152, 14, 40, 84, 68,252,  0, 66,121,127,223, 78,
		201, 34,225,240,124,191, 12, 42, 92,209,  1, 31,130,100, 28,224,
		161,249,110, 77, 87,144,181, 21, 86, 58,174,113,194,147,242, 50,
		 59, 48,250, 88,184,245, 45, 44,148, 73,154,230,149, 89,118,119,
		 79,229,239,117,189,179,254,155, 20,176,157,212, 36,123,234, 46,
		159,202,171, 67, 93, 62, 47,164,247, 15,137,235,216,160,200,133,
		140,172,192,221,131,111,218,210,153,219, 41, 72, 63,183, 39,  8,
		 98,168, 52,213,175,134,115, 90, 82, 16,226,220,251, 69,243,233,
		180,  5, 99, 60,  7,204,112,253,182,109, 25, 53, 94,187,190,178,
		 97, 49,126,  6, 64,143,169,173, 43,199, 33, 96, 76, 11,236, 23,
		166,104,141,246,125,217,122,238, 27, 22,228,146,165, 71,232,222
	};
	uint8_t h=seed;
	for(size_t i=0;i<nwords;i++) {
		h = Rand8[h ^ (uint8_t)(kmer[i]>>56)];
		h = Rand8[h ^ (uint8_t)(kmer[i]>>48) & 255];
		h = Rand8[h ^ (uint8_t)(kmer[i]>>40) & 255];
		h = Rand8[h ^ (uint8_t)(kmer[i]>>32) & 255];
		h = Rand8[h ^ (uint8_t)(kmer[i]>>24) & 255];
		h = Rand8[h ^ (uint8_t)(kmer[i]>>16) & 255];
		h = Rand8[h ^ (uint8_t)(kmer[i]>>8) & 255];
		h = Rand8[h ^ (uint8_t)kmer[i] & 255];
	}
	return h;
}

inline word_t* kmerizer::canonicalize(word_t *packed, word_t *rcpack) const {
	int cmp=0;
	for(size_t i=0;i<nwords;i++) {
		rcpack[i] = revcomp(packed[nwords-1-i]);
		if (rshift && i==nwords-1)
			rcpack[i] >>= rshift;
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

#endif
