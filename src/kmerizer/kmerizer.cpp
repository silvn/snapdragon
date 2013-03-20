#include "kmerizer.h"
#include <boost/thread.hpp>
#include <algorithm> // stl sort
#include <cstdlib> // qsort()
#include <cstdio>  // sprintf()
#include <cstring> // memcpy()
#include <sys/stat.h> // mkdir()

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
	state = READING;
	batches=0;
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
	serialize();
	
	if (batches>1) {
		int rc = mergeBatches();
		fprintf(stderr,"error merging batches\n");
		return 1;
	}
}

int kmerizer::histogram() {
	if (batches>1) {
		int rc = save();
		fprintf(stderr,"error saving\n");
		exit(1);
	}
	if (state == READING) {
		int rc = unique();
		if (rc != 0) {
			fprintf(stderr,"error uniqifying kmers\n");
			return 1;
		}
	}
	// merge the NBINS counts bvecs
	uint32_t todo = NBINS;
	size_t offset [NBINS];
	size_t done [NBINS];
	memset(offset,0,sizeof(uint32_t)*NBINS);
	memset(done,0,sizeof(uint32_t)*NBINS);
	uint32_t key=1; // min kmer frequency
	while (todo > 0) {
		uint32_t val = 0;
		for(size_t i=0;i<NBINS;i++) {
			if (done[i]==0) {
				if (kmer_freq[i][offset[i]] == key) {
					val += counts[i][offset[i]]->count();
					if (offset[i]>0) // undo the range encoding <=
						val -= counts[i][offset[i-1]]->count();
					offset[i]++;
					if(offset[i] == kmer_freq[i].size()) {
						done[i]=1;
						todo--;
					}
				}
			}
		}
		if (val>0) {
			printf("%u %u\n",key,val);
		}
		key++;
	}
	return 0;
}

void kmerizer::serialize() {
	// unique the batch
	int rc = unique();
	if (rc != 0) {
		fprintf(stderr,"error uniqifying a batch of kmers\n");
		exit(1);
	}

	// write to disk
	rc = writeBatch();
	if (rc != 0) {
		fprintf(stderr,"error writing a batch of kmers to disk\n");
		exit(1);
	}
	batches++;

	// zero the data
	memset(bin_tally,0,sizeof(uint32_t)*NBINS);
	for(size_t i=0;i<NBINS;i++) {
		delete counts[i]; // bvec32 destructor
	}
	state = READING;
}

int kmerizer::unique() {
	boost::thread_group tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j=(i+thread_bins>NBINS) ? NBINS : i+thread_bins;
		tg.create_thread(boost::bind(do_unique,i,j));
	}
	tg.join_all();
	state = QUERY;
}

int kmerizer::writeBatch() {
	boost::thread_group tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j=(i+thread_bins>NBINS) ? NBINS : i+thread_bins;
		tg.create_thread(boost::bind(do_writeBatch,i,j));
	}
	tg.join_all();
}

int kmerizer::mergeBatches() {
	boost::thread_group tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j=(i+thread_bins>NBINS) ? NBINS : i+thread_bins;
		tg.create_thread(boost::bind(do_mergeBatches,i,j));
	}
	tg.join_all();
}

void kmerizer::do_unique(const size_t from, const size_t to) {
	for(size_t bin=from; bin<to; bin++) {
		// sort
		qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers);
		// uniq
		uint32_t distinct = 0;
		vector<uint32_t> tally;
		tally.push_back(1); // first kmer
		for(size_t i=1;i<bin_tally[bin];i++) {
			word_t *ith = kmer_buf[bin] + i*nwords;
			if(compare_kmers(kmer_buf[bin] + distinct*nwords, ith) == 0)
				tally.back()++;
			else {
				distinct++;
				tally.push_back(1);
				memcpy(kmer_buf[bin] + distinct*nwords, ith, kmer_size);
			}
		}
		bin_tally[bin] = distinct+1;
		// create a bitmap index for the tally vector
		range_index(tally,kmer_freq[bin],counts[bin]);
	}
}

void kmerizer::do_writeBatch(const size_t from, const size_t to) {
	
}

void kmerizer::do_mergeBatches(const size_t from, const size_t to) {
	
}


// for each distinct value in the vec create a bitvector
// indexing the positions in the vec holding a value <= v
void range_index(vector<uint32_t> &vec, vector<uint32_t> &values, vector<bvec32*> &index) {
	// find the distinct values in vec
	values.clear();
	values.reserve(vec.size());
	copy(vec.begin(),vec.end(),values.begin());
	sort(values.begin(),values.end());
	vector<uint32_t>::iterator it;
	it = unique(values.begin(),values.end());
	values.resize(distance(values.begin(),it));
	// setup a vector for each range
	vector<uint32_t> range [values.size()];
	// iterate over the vec and push the offset onto each range
	for(size_t i=0;i<vec.size();i++) {
		it = find(values.begin(),values.end(),vec[i]);
		for(size_t j=0;j<=distance(values.begin(),it);j++) {
			range[j].push_back(i);
		}
	}
	// create a bitvector for each range
	index.resize(values.size());
	for(size_t i=0;i<range.size();i++) {
		index[i] = new bvec32(range[i]);
		fprintf(stderr,"index[%zi]->bytes()=%u\n",i,index[i]->bytes());
	}
}
