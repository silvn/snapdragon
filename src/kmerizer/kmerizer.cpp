#include "kmerizer.h"
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <bitset>
#include <algorithm> // stl sort
#include <cstdlib> // qsort()
#include <cstdio>  // sprintf()
#include <cstring> // memcpy()
#include <sys/stat.h> // mkdir()

// constructor
kmerizer::kmerizer(const size_t _k, const size_t _threads, char* _outdir, const char _mode) {
	k = _k;
	kmask = 0xFFFFFFFFFFFFFFFFULL;
	shiftlastby = 62;
	lshift=0;
	rshift=0;
	if (k%32 > 0) {
		kmask = (1ULL << (2*(k%32))) - 1;
		shiftlastby = 2*(k%32) - 2;
		lshift = 2*(k%32);
		rshift = 64 - lshift;
	}
	nwords = ((k-1)>>5)+1;
	kmer_size = nwords*sizeof(word_t);
	outdir = _outdir;
	mode = _mode;
	threads = _threads;
	thread_bins = NBINS/threads;
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

word_t* kmerizer::canonicalize(word_t *packed, word_t *rcpack) const {
	for(size_t i=0;i<nwords;i++) {
		rcpack[i] = packed[nwords-1-i];
	}
	if (lshift) {
		for(size_t i=0;i<nwords-1;i++) {
			rcpack[i] |= rcpack[i+1] << lshift;
			rcpack[i+1] >>= rshift;
		}
	}
	int cmp=0;
	for(size_t i=0;i<nwords;i++) {
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

void kmerizer::next_kmer(word_t* kmer, const char nucl) {
	// shift first kmer by 2 bits
	kmer[0] <<= 2;
	for(size_t w=1;w<nwords-1;w++) { // middle (full length) words
		kmer[w-1] |= kmer[w] >> 62; // move the top two bits
		kmer[w] <<= 2; // make room for the next word
	}
	// last word (if not the first)
	if(nwords>1) {
		kmer[nwords-2] |= (kmer[nwords-1] >> shiftlastby);
		kmer[nwords-1] <<= 2;
	}
	kmer[nwords-1] |= twobit(nucl);
	kmer[nwords-1] &= kmask;
}

void kmerizer::addSequence(const char* seq, const int length) {
	if (length < k) return;
	word_t packed [nwords];
	word_t rcpack [nwords];
	memset(packed,0,kmer_size);

	for(size_t i=0;i<k;i++)
		next_kmer(packed,seq[i]);

	word_t *kmer = packed;
	size_t bin;
	if (mode == CANONICAL)
		kmer = canonicalize(packed,rcpack);
	if (mode == BOTH) {
		kmer = canonicalize(packed,rcpack);
		bin = hashkmer(packed,0);
		memcpy(kmer_buf[bin] + nwords*bin_tally[bin], packed, kmer_size);
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
		next_kmer(packed, seq[i]);
		if (mode == CANONICAL)
			kmer = canonicalize(packed,rcpack);
		if (mode == BOTH) {
			kmer = canonicalize(packed,rcpack);
			bin = hashkmer(packed,0);
			memcpy(kmer_buf[bin] + nwords*bin_tally[bin], packed, kmer_size);
			bin_tally[bin]++;
			if (bin_tally[bin] == max_kmers_per_bin) serialize();
			kmer = rcpack;
		}
		bin = hashkmer(kmer,0);
		memcpy(kmer_buf[bin] + nwords*bin_tally[bin], kmer, kmer_size);
		bin_tally[bin]++;
		if (bin_tally[bin] == max_kmers_per_bin) serialize();
	}
}

void kmerizer::save() {
	serialize();
	
	if (batches>1)
		mergeBatches();
	else
		for(size_t bin=0;bin<NBINS;bin++) {
			char ofname [100];
			char nfname [100];
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

void kmerizer::load() {
	boost::thread_group tg;
	for(size_t i = 0; i < NBINS; i += thread_bins) {
		size_t j = (i + thread_bins > NBINS) ? NBINS : i + thread_bins;
		tg.create_thread(boost::bind(&kmerizer::do_loadIndex, this, i, j));
	}
	tg.join_all();
	state = QUERY;
}

// given one kmer, pack it, canonicalize it, hash it, find it
uint32_t kmerizer::find(char* seq) {
	// pack it
	word_t packed [nwords];
	word_t rcpack [nwords];
	memset(packed,0,kmer_size);
	for(size_t i=0;i<k;i++)
		next_kmer(packed,seq[i]);
	word_t *kmer = packed;

	// canonicalize it
	if (mode == CANONICAL) {
		kmer = canonicalize(packed,rcpack);
	}

	//hash it
	size_t bin = hashkmer(kmer,0);

	// find it
	return find(kmer,bin);
}

uint32_t kmerizer::find(word_t *kmer, size_t bin) {
	if (kmer_buf[bin]==NULL) {
		// try to mmap the file into the kmer_buf
		// if mmap fails, do a linear scan of the file and return the associated count
	}
	// we've got the data in kmer_buf[bin]
}

void kmerizer::histogram() {
	if (batches>1)
		save();

	if (state == READING)
		uniqify();

	// merge the NBINS counts bvecs
	uint32_t todo = NBINS;
	size_t offset [NBINS];
	size_t done [NBINS];
	memset(offset,0,sizeof(size_t)*NBINS);
	memset(done,0,sizeof(size_t)*NBINS);
	uint32_t key=1; // min kmer frequency
	while (todo > 0) {
		uint32_t val = 0;
		for (size_t i=0;i<NBINS;i++) {
			if (done[i]==0) {
				if (kmer_freq[i][offset[i]] == key) {
					val += counts[i][offset[i]]->cnt();
					if (offset[i]<kmer_freq[i].size()-1) // undo the range encoding <=
						val -= counts[i][offset[i]+1]->cnt();
					offset[i]++;
					if (offset[i] == kmer_freq[i].size()) {
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

void kmerizer::serialize() {
	batches++;
	// unique the batch
	uniqify();
	// write to disk
	writeBatch();

	// zero the data
	memset(bin_tally,0,sizeof(uint32_t)*NBINS);
	for(size_t i=0;i<NBINS;i++) {
		for(size_t j=0; j<counts[i].size(); j++)
			delete counts[i][j]; // bvec32 destructor
	}
	state = READING;
}

void kmerizer::uniqify() {
	boost::thread_group tg;
	for(size_t i = 0; i < NBINS; i += thread_bins) {
		size_t j = (i + thread_bins > NBINS) ? NBINS : i + thread_bins;
		tg.create_thread(boost::bind(&kmerizer::do_unique, this, i, j));
	}
	tg.join_all();
	state = QUERY;
}

void kmerizer::writeBatch() {
	boost::thread_group tg;
	for(size_t i = 0; i < NBINS; i += thread_bins) {
		size_t j = (i + thread_bins > NBINS) ? NBINS : i + thread_bins;
		tg.create_thread(boost::bind(&kmerizer::do_writeBatch, this, i, j));
	}
	tg.join_all();
}

void kmerizer::mergeBatches() {
	return;
	boost::thread_group tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j=(i+thread_bins>NBINS) ? NBINS : i+thread_bins;
		tg.create_thread(boost::bind(&kmerizer::do_mergeBatches, this, i, j));
	}
	tg.join_all();
	batches=1;
}


int compare_kmers1(const void *k1, const void *k2) {
	if (*(word_t*)k1 < *(word_t*)k2) return -1;
	if (*(word_t*)k1 > *(word_t*)k2) return 1;
	return 0;
}
int compare_kmers2(const void *k1, const void *k2) {
	for(size_t i=0;i<2;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}
int compare_kmers3(const void *k1, const void *k2) {
	for(size_t i=0;i<3;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}
int compare_kmers4(const void *k1, const void *k2) {
	for(size_t i=0;i<4;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}
int compare_kmers5(const void *k1, const void *k2) {
	for(size_t i=0;i<5;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}
int compare_kmers6(const void *k1, const void *k2) {
	for(size_t i=0;i<6;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}
int compare_kmers7(const void *k1, const void *k2) {
	for(size_t i=0;i<7;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}
int compare_kmers8(const void *k1, const void *k2) {
	for(size_t i=0;i<8;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}

int kmercmp(const void *k1, const void *k2, size_t nwords) {
	for(size_t i=0;i<nwords;i++) {
		if (*((word_t*)k1+i) < *((word_t*)k2+i)) return -1;
		if (*((word_t*)k1+i) > *((word_t*)k2+i)) return 1;
	}
	return 0;
}

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
			if(kmercmp(kmer_buf[bin] + distinct*nwords, ith, nwords) == 0)
				tally.back()++;
			else {
				distinct++;
				tally.push_back(1);
				memcpy(kmer_buf[bin] + distinct*nwords, ith, kmer_size);
			}
		}
//		fprintf(stderr,"uniqify[%zi] reduced from %u - %u = %u\n",bin,bin_tally[bin],distinct+1,bin_tally[bin] - distinct - 1);
		bin_tally[bin] = distinct+1;
		// create a bitmap index for the tally vector
		range_index(tally,kmer_freq[bin],counts[bin]);
	}
}

void kmerizer::print_kmer(word_t *kmer) {
	size_t bpw = 8*sizeof(word_t);
	for(size_t j=0;j<nwords;j++) {
		fprintf(stderr," ");
		for(size_t b=0;b<bpw;b++) {
			fprintf(stderr,"%d",((*(kmer+j) & (1ULL << (bpw-b-1)))) ? 1 : 0);
		}
	}
	fprintf(stderr,"\n");
}
void kmerizer::bit_slice(word_t *kmers, const size_t n, bvec32 **kmer_slices, size_t nbits) {
	for(size_t i=0;i<nbits;i++) {
		kmer_slices[i] = new bvec32();
	}
	size_t bword [nbits]; // which word does the bth bit live in?
	word_t bpos [nbits]; // what's the position of the bth bit?
	word_t brun [nbits]; // how many same bits have we seen in the bth bit?
	bool bbit [nbits]; // which bit did we have here last time?
	size_t bpw = 8*sizeof(word_t); // bits per word
	for(size_t b=0;b<nbits;b++) {
		bword[b] = b/bpw;
		bpos[b] = 1ULL << (bpw - b - 1);
		brun[b] = 0;
		bbit[b] = false;
	}
	for(size_t i=0;i<n;i++) {
		for(size_t b=0;b<nbits;b++) {
			// check the bth bit
			if (kmers[i*nwords + bword[b]] & bpos[b]) { // bth bit is set
				// extend the 1-run or appendFill?
				if (bbit[b])
					brun[b]++;
				else {
					if(brun[b]>0) 
						kmer_slices[b]->appendFill(bbit[b],brun[b]);
					brun[b] = 0;
					bbit[b]=true;
				}
			}
			else { // bth bit is not set
				// extend the 0-run or appendFill?
				if (bbit[b]) {
					if(brun[b]>0)
						kmer_slices[b]->appendFill(bbit[b],brun[b]);
					brun[b]=0;
					bbit[b]=false;
				}
				else
					brun[b]++;
			}
		}
	}
	// append the last runs
	for(size_t b=0; b<nbits; b++)
		if(brun[b]>0)
			kmer_slices[b]->appendFill(bbit[b],brun[b]);
}

void kmerizer::do_writeBatch(const size_t from, const size_t to) {
	for(size_t bin=from; bin<to; bin++) {
		// convert kmer_buf[bin] into a bit-sliced bitmap index
//		fprintf(stderr,"bin_tally[%zi]: %u\n",bin,bin_tally[bin]);
		const size_t nbits = 8*kmer_size;
		bvec32* kmer_slices [nbits];
		bit_slice(kmer_buf[bin],bin_tally[bin],kmer_slices,nbits);
		
		FILE *fp;
		if (0) {
			// open output file for kmer_buf
			char kmer_file [100];
			sprintf(kmer_file,"%s/%zi-mers.%zi.%zi.uncompressed",outdir,k,bin,batches);
			fp = fopen(kmer_file, "wb");
			fwrite(kmer_buf[bin],kmer_size,bin_tally[bin],fp);
			fclose(fp);
		}
		if (1) {
			// write the kmer_slices to a file
			char kmer_file [100];
			sprintf(kmer_file,"%s/%zi-mers.%zi.%zi",outdir,k,bin,batches);
			fp = fopen(kmer_file, "wb");
			size_t n_slices = 8*sizeof(word_t);
			fwrite(&n_slices,sizeof(size_t),1,fp);
			// fprintf(stderr,"bin_tally[%zi]: %u\n",bin,bin_tally[bin]);
			for(size_t b=0;b<n_slices;b++) {
				uint32_t *buf;
				size_t bytes = kmer_slices[b]->dump(&buf);
				// fprintf(stderr,"%zi: %zi,",b,bytes);
				fwrite(&bytes,sizeof(size_t),1,fp);
				fwrite(buf,1,bytes,fp);
				free(buf);
			}
//			fprintf(stderr,"\n");
			fclose(fp);
		}
		// open output file for counts
		char counts_file [100];
		sprintf(counts_file,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,batches);
		fp = fopen(counts_file, "wb");
		// first write the number of distinct values
		// then write the distinct values
		// for each distinct value, write out the bvec (size,count,rle,words.size(),words)
		size_t n_distinct = kmer_freq[bin].size();
		fwrite(&n_distinct,sizeof(size_t),1,fp);
		fwrite(kmer_freq[bin].data(),sizeof(uint32_t),n_distinct,fp);
		for(size_t i=0;i<n_distinct;i++) {
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

void kmerizer::read_index(const char* idxfile, vector<uint32_t> &values, vector<bvec32*> &index) {
	// open idxfile
	FILE *fp;
	fp = fopen(idxfile,"rb");
	// fread to repopulate values and bvecs
	size_t n_distinct;
	fread(&n_distinct,sizeof(size_t),1,fp);
	values.resize(n_distinct);
	fread(values.data(),sizeof(uint32_t),n_distinct,fp);
	index.resize(n_distinct);
	for(size_t i=0;i<n_distinct;i++) {
		size_t bytes;
		fread(&bytes,sizeof(size_t),1,fp);
		size_t words = bytes/sizeof(uint32_t);
		uint32_t buf [words];
		fread(buf,1,bytes,fp);
		index[i] = new bvec32(buf);
	}
	fclose(fp);
}

void kmerizer::do_loadIndex(const size_t from, const size_t to) {
	for(size_t bin=from; bin<to;bin++) {
		char fname [100];
		sprintf(fname,"%s/%zi-mers.%zi.idx",outdir,k,bin);
		read_index(fname,kmer_freq[bin],counts[bin]);
	}
}

void kmerizer::do_mergeBatches(const size_t from, const size_t to) {
	for(size_t bin=from; bin<to; bin++) {
		// read the counts for each batch
		// and open the files of sorted distinct kmers
		vector<bvec32*> batch_counts [batches];
		vector<uint32_t> batch_values [batches];
		FILE *kmer_fp [batches];
		for(size_t i=0;i<batches;i++) {
			char fname [100];
			sprintf(fname,"%s/%zi-mers.%zi.%zi",outdir,k,bin,i+1);
			kmer_fp[i] = fopen(fname, "rb");
			sprintf(fname,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,i+1);
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
			char distinct_seq [k];
			unpack(distinct,distinct_seq);
//			fprintf(stderr,"mindex: %zi distinct: %s\n",mindex,distinct_seq);
			size_t rc = fread(kmers + mindex*nwords,kmer_size,1,kmer_fp[mindex]);
			if (rc != 1) {
				todo--;
				btally[mindex]=0;
				if (todo==0) break;
			}
			else
				btally[mindex] = pos2value(offset[mindex]++,batch_values[mindex],batch_counts[mindex]);
			mindex = find_min(kmers,btally);
			if (kmercmp(kmers + mindex*nwords,distinct,nwords)==0)
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
			size_t bytes = counts[bin][i]->dump(&buf);
			fwrite(&bytes,sizeof(size_t),1,ofp);
			fwrite(counts[bin][i]->get_words().data(),1,counts[bin][i]->bytes(),ofp);
		}
		fclose(ofp);
		// delete input files
		for(size_t i=0;i<batches;i++) {
			sprintf(fname,"%s/%zi-mers.%zi.%zi",outdir,k,bin,i+1);
			if (remove(fname) != 0) perror("error deleting file");
			sprintf(fname,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,i+1);
			if (remove(fname) != 0) perror("error deleting file");
		}
	}
}

inline size_t kmerizer::find_min(const word_t* kmers, const uint32_t* kcounts) {
	size_t mindex = 0;
	while (kcounts[mindex] == 0) mindex++;
	for(size_t i=mindex+1;i<batches;i++)
		if (kcounts[mindex] != 0)
			if (kmercmp(kmers + i*nwords,kmers + mindex*nwords,nwords) < 0)
				mindex = i;

	return mindex;
}

// for each distinct value in the vec create a bitvector
// indexing the positions in the vec holding a value <= v
void kmerizer::range_index(vector<uint32_t> &vec, vector<uint32_t> &values, vector<bvec32*> &index) {
	// find the distinct values in vec
	values.clear();
	// use a bitset to mark the distinct values (from 0-255)
	// and another vector for overflow (which we'll sort/uniq later)
	bitset<256> mybits;
	vector<uint32_t> overflow;
	vector<uint32_t>::iterator it;
	for(it = vec.begin(); it < vec.end(); ++it) {
		if (*it >= 256) 
			overflow.push_back(*it);
		else
			mybits[*it]=1;
	}
	// populate values vector with set bits in mybits
	for(uint32_t i=0;i<256;i++)
		if (mybits.test(i))
			values.push_back(i);
	
	if (overflow.size() > 0) {
		sort(overflow.begin(),overflow.end());
		it = unique(overflow.begin(),overflow.end());
		values.insert(values.end(),overflow.begin(),it);
	}

	// setup a vector for each range
	vector<uint32_t> vrange [values.size()];
	// iterate over the vec and push the offset onto each range
	for(size_t i=0;i<vec.size();i++) {
		it = lower_bound(values.begin(),values.end(),vec[i]);
		for(size_t j=0;j<=(it - values.begin());j++)
			vrange[j].push_back(i);
	}
	// create a bitvector for each range
	index.resize(values.size());
	for(size_t i=0;i<values.size();i++)
		index[i] = new bvec32(vrange[i]);
}


uint32_t kmerizer::pos2value(size_t pos, vector<uint32_t> &values, vector<bvec32*> &index) {
	// lookup the value in the pos bit
	// find the first bvec where this bit is set
	for(size_t i=0;i<values.size();i++)
		if (index[i]->find(pos))
			return values[i];
	return 0;
}
