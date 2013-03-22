#include "kmerizer.h"
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <algorithm> // stl sort
#include <cstdlib> // qsort()
#include <cstdio>  // sprintf()
#include <cstring> // memcpy()
#include <sys/stat.h> // mkdir()

// constructor
kmerizer::kmerizer(const size_t _k, const size_t _threads, char* _outdir, const char _mode) {
	k = _k;
	kmask = 0xFFFFFFFFFFFFFFFF;
	shiftlastby = 62;
	if (k<32) {
		kmask = (1ULL << (k*2)) - 1;
		rshift = 64 - k*2;
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
	thread_bins = threads/NBINS;
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
		int rc = uniqify();
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
					val += counts[i][offset[i]]->cnt();
					if (offset[i]>0) // undo the range encoding <=
						val -= counts[i][offset[i-1]]->cnt();
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
	batches++;
	// unique the batch
	int rc = uniqify();
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

	// zero the data
	memset(bin_tally,0,sizeof(uint32_t)*NBINS);
	for(size_t i=0;i<NBINS;i++) {
		for(size_t j=0; j<counts[i].size(); j++)
			delete counts[i][j]; // bvec32 destructor
	}
	state = READING;
}

int kmerizer::uniqify() {
	boost::thread_group tg;
	for(size_t i = 0; i < NBINS; i += thread_bins) {
		size_t j = (i + thread_bins > NBINS) ? NBINS : i + thread_bins;
		tg.create_thread(boost::bind(&kmerizer::do_unique, this, i, j));
	}
	tg.join_all();
	state = QUERY;
}

int kmerizer::writeBatch() {
	boost::thread_group tg;
	for(size_t i = 0; i < NBINS; i += thread_bins) {
		size_t j = (i + thread_bins > NBINS) ? NBINS : i + thread_bins;
		tg.create_thread(boost::bind(&kmerizer::do_writeBatch, this, i, j));
	}
	tg.join_all();
}

int kmerizer::mergeBatches() {
	boost::thread_group tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j=(i+thread_bins>NBINS) ? NBINS : i+thread_bins;
		tg.create_thread(boost::bind(&kmerizer::do_mergeBatches, this, i, j));
	}
	tg.join_all();
	batches=1;
}


int compare_kmers1(const void *k1, const void *k2) { return memcmp(k1,k2,sizeof(word_t)); }
int compare_kmers2(const void *k1, const void *k2) { return memcmp(k1,k2,2*sizeof(word_t)); }
int compare_kmers3(const void *k1, const void *k2) { return memcmp(k1,k2,3*sizeof(word_t)); }
int compare_kmers4(const void *k1, const void *k2) { return memcmp(k1,k2,4*sizeof(word_t)); }
int compare_kmers5(const void *k1, const void *k2) { return memcmp(k1,k2,5*sizeof(word_t)); }
int compare_kmers6(const void *k1, const void *k2) { return memcmp(k1,k2,6*sizeof(word_t)); }
int compare_kmers7(const void *k1, const void *k2) { return memcmp(k1,k2,7*sizeof(word_t)); }
int compare_kmers8(const void *k1, const void *k2) { return memcmp(k1,k2,8*sizeof(word_t)); }

void kmerizer::do_unique(const size_t from, const size_t to) {
	for(size_t bin=from; bin<to; bin++) {
		// sort
		switch(nwords) {
			case 1:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers1);
			case 2:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers2);
			case 3:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers3);
			case 4:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers4);
			case 5:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers5);
			case 6:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers6);
			case 7:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers7);
			case 8:
				qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers8);
		}
		// uniq
		uint32_t distinct = 0;
		vector<uint32_t> tally;
		tally.push_back(1); // first kmer
		for(size_t i=1;i<bin_tally[bin];i++) {
			word_t *ith = kmer_buf[bin] + i*nwords;
			if(memcmp(kmer_buf[bin] + distinct*nwords, ith, kmer_size) == 0)
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
	for(size_t bin=from; bin<to; bin++) {
		// open output file for kmer_buf
		char kmer_file [100];
		sprintf(kmer_file,"%s/%zi-mers.%zi.%zi",outdir,k,bin,batches);
		FILE *fp;
		fp = fopen(kmer_file, "wb");
		fwrite(kmer_buf[bin],kmer_size,bin_tally[bin],fp);
		fclose(fp);

		// open output file for counts
		char counts_file [100];
		sprintf(counts_file,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,batches);
		fp = fopen(counts_file, "wb");
		// first write the number of distinct values
		// then write the distinct values
		// for each distinct value, write out the bvec (size,count,rle,words.size(),words)
		int n_distinct = kmer_freq[bin].size();
		fwrite(&n_distinct,sizeof(int),1,fp);
		fwrite(kmer_freq[bin].data(),sizeof(uint32_t),n_distinct,fp);
		for(size_t i=0;i<n_distinct;i++) {
			uint32_t* buf;
			size_t bytes = counts[bin][i]->dump(buf);
			// how many bytes are we writing
			fwrite(&bytes,sizeof(size_t),1,fp);
			// write them
			fwrite(buf,1,bytes,fp);
		}
		fclose(fp);
	}
}

void kmerizer::read_index(const char* idxfile, vector<uint32_t> &values, vector<bvec32*> &index) {
	// open idxfile
	FILE *fp;
	fp = fopen(idxfile,"rb");
	// fread to repopulate values and bvecs
	int n_distinct;
	fread(&n_distinct,sizeof(int),1,fp);
	values.resize(n_distinct);
	fread(values.data(),sizeof(uint32_t),n_distinct,fp);
	for(size_t i=0;i<n_distinct;i++) {
		size_t bytes;
		fread(&bytes,sizeof(size_t),1,fp);
		size_t words = bytes/sizeof(uint32_t);
		uint32_t buf [words];
		fread(buf,sizeof(uint32_t),words,fp);
		index[i] = new bvec32(buf);
	}
}

void kmerizer::do_mergeBatches(const size_t from, const size_t to) {
	for(size_t bin=from; bin<to; bin++) {
		// read the counts for each batch
		// and open the files of sorted distinct kmers
		vector<bvec32*> batch_counts [batches];
		vector<uint32_t> batch_values [batches];
		FILE *kmer_fp [batches];
		for(size_t i=1;i<=batches;i++) {
			char fname [100];
			sprintf(fname,"%s/%zi-mers.%zi.%zi",outdir,k,bin,i);
			kmer_fp[i-1] = fopen(fname, "rb");
			sprintf(fname,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,i);
			read_index(fname,batch_values[i],batch_counts[i]);
		}
		// open an output file for the merged distinct kmers
		FILE *ofp;
		char fname [100];
		sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);
		ofp = fopen(fname,"wb");
		vector<uint32_t> tally;
		// read the first kmer from each batch
		word_t kmers [batches*nwords];
		uint32_t btally [batches];
		uint32_t todo = batches;
		size_t offset [batches];
		memset(offset,0,sizeof(size_t)*batches);
		for(size_t i=0;i<batches;i++) {
			size_t rc = fread(kmers + i*nwords, kmer_size,1,kmer_fp[i]);
			if (rc != 1) {
				todo--;
				btally[i]=0;
			}
			else
				btally[i] = pos2value(offset[i]++,batch_values[i],batch_counts[i]);
		}
		// choose min
		size_t mindex = find_min(kmers,btally);
		word_t distinct [nwords];
		memcpy(distinct,kmers + mindex*nwords, kmer_size);
		tally.push_back(btally[mindex]);
		// replace min
		// iterate until there's nothing left to do
		while(todo>0) {
			size_t rc = fread(kmers + mindex*nwords,kmer_size,1,kmer_fp[mindex]);
			if (rc != 1) {
				todo--;
				btally[mindex]=0;
				if (todo==0) break;
			}
			else
				btally[mindex] = pos2value(offset[mindex]++,batch_values[mindex],batch_counts[mindex]);
			mindex = find_min(kmers,btally);
			if (memcmp(kmers + mindex*nwords,distinct,kmer_size)==0)
				tally.back() += btally[mindex];
			else {
				size_t rc = fwrite(distinct,kmer_size,1,ofp);
				if (rc != 1) {
					fprintf(stderr,"error writing distinct kmer to output file\n");
					exit(1);
				}
				memcpy(distinct,kmers+mindex*nwords,kmer_size);
				tally.push_back(btally[mindex]);
			}
		}
		range_index(tally,kmer_freq[bin],counts[bin]);
		// close output file
		fclose(ofp);
		// close input files
		for(size_t i=0;i<batches;i++)
			fclose(kmer_fp[i]);
		
		char counts_file [100];
		sprintf(counts_file,"%s/%zi-mers.%zi.idx",outdir,k,bin);
		ofp = fopen(counts_file, "wb");
		// first write the number of distinct values
		// then write the distinct values
		// for each distinct value, write out the bvec (size,count,rle,words.size(),words)
		int n_distinct = kmer_freq[bin].size();
		fwrite(&n_distinct,sizeof(int),1,ofp);
		fwrite(kmer_freq[bin].data(),sizeof(uint32_t),n_distinct,ofp);
		for(size_t i=0;i<n_distinct;i++) {
			// I want this DIY serialization/deserialization in bvec
			uint32_t *buf;
			size_t bytes = counts[bin][i]->dump(buf);
			fwrite(&bytes,sizeof(size_t),1,ofp);
			fwrite(counts[bin][i]->get_words().data(),1,counts[bin][i]->bytes(),ofp);
		}
		fclose(ofp);
		// delete input files
		
	}
}

inline size_t kmerizer::find_min(const word_t* kmers, const uint32_t* kcounts) {
	size_t mindex = 0;
	while (kcounts[mindex] == 0) mindex++;
	for(size_t i=mindex+1;i<batches;i++)
		if (kcounts[mindex] != 0)
			if (memcmp(kmers + i*nwords,kmers + mindex*nwords,kmer_size) < 0)
				mindex = i;

	return mindex;
}

// for each distinct value in the vec create a bitvector
// indexing the positions in the vec holding a value <= v
void kmerizer::range_index(vector<uint32_t> &vec, vector<uint32_t> &values, vector<bvec32*> &index) {
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
	for(size_t i=0;i<values.size();i++) {
		index[i] = new bvec32(range[i]);
		fprintf(stderr,"index[%zi]->bytes()=%u\n",i,index[i]->bytes());
	}
}


uint32_t
kmerizer::pos2value(
	size_t pos,
	vector<uint32_t> &values,
	vector<bvec32*> &index
) {
	// lookup the value in the pos bit
	// find the first bvec where this bit is set
	return 0;
}

int main (int argc, char const *argv[]) {
	return 0;
}