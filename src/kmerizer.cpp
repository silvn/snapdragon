#include "kmerizer.h"
#include <zlib.h>
#include "kseq.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <boost/thread.hpp>

// constructor
kmerizer::kmerizer(const size_t _k, const size_t _threads, const char* _outdir, const char _mode) {
	k = _k;
	kmask = 0xFFFFFFFFFFFFFFFF;
	shiftlastby = 62;
	if (k<32) {
		kmask = (1ULL << (k*2)) - 1;
		rshift = 64 - k*2
	}
	if (k>32 & k%32 > 0) {
		kmask = (1ULL << (2*(k%32))) - 1;
		shiftlastby = 2*(k%32) - 2;
		rshift = 64 - 2*(k%32);
	}
	nwords = ((k-1)>>5)+1;
	kmer_size = nwords*sizeof(word_t);
	outdir = _outdir;
	mode = _mode;
	threads = _threads;
	thread_bins = threads/NBINS:
}

int kmerizer::allocate(const size_t maximem) {
	memset(bin_tally,0,sizeof(uint32_t)*NBINS);
	max_kmers_per_bin = maximem/kmer_size/NBINS;
	for(size_t i=0; i<NBINS; i++) {
		kmer_buf[i] = (word_t*) calloc(max_kmers_per_bin, kmer_size);
		if (kmer_buf[i]==NULL) return 1;
	}
	return 0;
}

int kmerizer::addSequence(const char* seq, const int length) {
	if (length < k) return 1;
	word_t packed [nwords];
	word_t rcpack [nwords];
	memset(packed,0,kmer_size);
	packed[0] = twobit(seq[0]);
	for(size_t i=1;i<k;i++) {
		size_t word = i>>5;
		packed[word] <<= 2;
		packed[word] |= twobit(seq[i]);
	}
	word_t *kmer = packed;
	size_t bin;
	if (mode == CANONICAL) {
		kmer = canonicalize(packed,rcpack);
	}
	if (mode == BOTH) {
		kmer = canonicalize(packed,rcpack);
		bin = hashkmer(packed,0);
		memcpy(kmer_buf[bin] + nwords*bin_tally[bin], kmer, kmer_size);
		bin_tally[bin]++;
		if (bin_tally[bin] == max_kmers_per_bin) serialize();
		kmer = rcpack;
	}
	bin = hashkmer(kmer,0);
	memcpy(kmer_buf[bin] + nwords*bin_tally[bin], kmer, kmer_size);
	bin_tally[bin]++;
	if (bin_tally[bin] == max_kmers_per_bin) serialize();
	// pack the rest of the sequence
	for(size_t i=k; i<(size_t)length;i++) {
		packed[0] <<= 2;
		if (k<32) packed[0] &= kmask;
		for(size_t word=1;word<nwords-1;word++) {
			packed[word-1] |= packed[word] >> 62;
			packed[word] <<= 2;
		}
		if (nwords>1) {
			packed[nwords-2] |= packed[nwords-1] >> shiftlastby;
			packed[nwords-1] &= kmask;
		}
		packed[nwords-1] |= twobit(seq[i]);
		if (mode == CANONICAL) {
			kmer = canonicalize(packed,rcpack);
		}
		if (mode == BOTH) {
			kmer = canonicalize(packed,rcpack);
			bin = hashkmer(packed,0);
			memcpy(kmer_buf[bin] + nwords*bin_tally[bin], kmer, kmer_size);
			bin_tally[bin]++;
			if (bin_tally[bin] == max_kmers_per_bin) serialize();
			kmer = rcpack;
		}
		bin = hashkmer(kmer,0);
		memcpy(kmer_buf[bin] + nwords*bin_tally[bin], kmer, kmer_size);
		bin_tally[bin]++;
		if (bin_tally[bin] == max_kmers_per_bin) serialize();
	}
	return 0;
}

int kmerizer::save() {
	
}

int kmerizer::histogram() {
	
}

void serialize() {
	// unique the batch
	// write to disk
	// zero the data
	memset(bin_tally,0,sizeof(uint32_t)*NBINS);
}

int unique() {
	boost::thread_group tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j=(i+thread_bins>NBINS) ? NBINS : i+thread_bins;
		tg.create_thread(boost::bind(sortuniq,i,j));
	}
	tg.join_all();
}

void sortuniq(const size_t from, const size_t to) {
	
}