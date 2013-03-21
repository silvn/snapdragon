#include <zlib.h>
#include <vector>
#include <algorithm> // stl sort
#include <stdlib.h> // qsort()
#include <stdio.h>  // sprintf()
#include <string.h> // memcpy()
#include <sys/stat.h> // mkdir()
#include "kseq.h"
#include <boost/thread.hpp>
using namespace std;

#define NBINS 256

size_t nwords;
size_t kmer_size;
size_t threads;
size_t thread_bins;
char* outprefix;

inline uint8_t hashkmer(uint64_t *kmer, uint8_t seed) {
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

inline uint64_t twobit(uint64_t val) {
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
	return (uint64_t)table[val];
}

inline uint64_t revcomp(uint64_t val) {
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
		((uint64_t)table[val&0xFFUL]<<56) |
		((uint64_t)table[(val>>8)&0xFFUL]<<48) |
		((uint64_t)table[(val>>16)&0xFFUL]<<40) |
		((uint64_t)table[(val>>24)&0xFFUL]<<32) |
		((uint64_t)table[(val>>32)&0xFFUL]<<24) |
		((uint64_t)table[(val>>40)&0xFFUL]<<16) |
		((uint64_t)table[(val>>48)&0xFFUL]<<8) |
		((uint64_t)table[(val>>56)&0xFFUL]);
}

int compare_kmers (const void* v1, const void *v2) {
	return memcmp(v1,v2,kmer_size); // system developers probably did a better job...
	uint64_t *k1 = (uint64_t*)v1;
	uint64_t *k2 = (uint64_t*)v2;
	for(size_t i=0;i<nwords;i++) {
		if (k1[i] > k2[i])
			return 1;
		else if (k1[i] < k2[i])
			return -1;
	}
	return 0;
}

void sortuniq(const size_t from, const size_t to, vector<uint64_t*> &kmer_buf, vector<size_t> &bin_tally, const size_t batch, const size_t k) {
	for (size_t bin=from;bin<to;bin++) {
		qsort(kmer_buf[bin], bin_tally[bin], kmer_size, compare_kmers);
		// open output files for this batch
		gzFile kmer_fp;
		char kmer_file [100];
		sprintf(kmer_file,"%s/distinct_%zi-mers.%zi.%zi",outprefix,k,bin,batch);
		kmer_fp = gzopen(kmer_file, "wb0");
		gzbuffer(kmer_fp,65536);
		gzFile kcount_fp;
		char kcount_file [100];
		sprintf(kcount_file,"%s/distinct_%zi-mer_counts.%zi.%zi",outprefix,k,bin,batch);
		kcount_fp = gzopen(kcount_file, "wb");
		gzbuffer(kcount_fp,65536);

		// we want to write gzipped uniq kmers and counts
		uint32_t distinct [nwords];
		memcpy(distinct,kmer_buf[bin],kmer_size);
		uint32_t count = 1;
		for(size_t i=1;i<bin_tally[bin];i++) {
			// is ith kmer same as distinct?
			uint64_t *ith = kmer_buf[bin] + i*nwords;
			if (compare_kmers(distinct,ith) == 0) // same kmer
				count++;
			else {
				gzwrite(kmer_fp, distinct, kmer_size);
				gzwrite(kcount_fp, &count, sizeof(uint32_t));
				memcpy(distinct,ith,kmer_size);
				count=1;
			}
		}
		// write last distinct kmer and count
		gzwrite(kmer_fp, distinct, kmer_size);
		gzwrite(kcount_fp, &count, sizeof(uint32_t));
		// close output files
		gzclose(kmer_fp);
		gzclose(kcount_fp);
	}
}

inline size_t find_min(const uint64_t* kmers,const uint32_t* kcounts, size_t batches) {
	size_t mindex = 0;
	while (kcounts[mindex] == 0) mindex++;
	for(size_t i=mindex+1;i<batches;i++)
		if (kcounts[mindex] != 0) {
			const uint64_t *k1 = kmers + i*nwords;
			const uint64_t *k2 = kmers + mindex*nwords;
			for(size_t j=0;j<nwords;j++) {
				if(k1[j] > k2[j])
					break;
				else if(k1[j] < k2[j]) {
					mindex = i;
					break;
				}
			}
		}
	return mindex;
}

void mergecounts(const size_t from, const size_t to, const size_t batches, size_t k) {
	fprintf(stderr,"mergecounts(%zi,%zi,%zi)\n",from,to,batches);
	// input files
	for(size_t bin=from;bin<to;bin++) {
		gzFile kmer_fp [batches];
		gzFile kcount_fp [batches];
		for(size_t i=0;i<batches;i++) {
			char kmer_file [100];
			sprintf(kmer_file,"%s/distinct_%zi-mers.%zi.%zi",outprefix,k,bin,i+1);
			kmer_fp[i] = gzopen(kmer_file, "rb");
			gzbuffer(kmer_fp[i],65536);
			char kcount_file [100];
			sprintf(kcount_file,"%s/distinct_%zi-mer_counts.%zi.%zi",outprefix,k,bin,i+1);
			kcount_fp[i] = gzopen(kcount_file, "rb");
			gzbuffer(kcount_fp[i],65536);
		}
		// open a pair of output files
		gzFile kmer_ofp;
		char kmer_ofile [100];
		sprintf(kmer_ofile,"%s/distinct_%zi-mers.%zi.0",outprefix,k,bin);
		kmer_ofp = gzopen(kmer_ofile, "wb0");
		gzbuffer(kmer_ofile,65536);
		gzFile kcount_ofp;
		char kcount_ofile [100];
		sprintf(kcount_ofile,"%s/distinct_%zi-mer_counts.%zi.0",outprefix,k,bin);
		kcount_ofp = gzopen(kcount_ofile, "wb");
		gzbuffer(kcount_ofp,65536);

		// read one kmer from each batch
		uint64_t kmers [batches*nwords];
		uint32_t kcounts [batches];
		uint32_t todo = batches;
		for(size_t i=0;i<batches;i++) {
			int rc = gzread(kmer_fp[i],kmers+i*nwords,kmer_size);
			if (rc == 0) {
				todo--;
				kcounts[i]=0;
			}
			else if (rc < 0) {
				fprintf(stderr,"some error code from gzread\n");
				todo--;
				kcounts[i]=0;
			}
			else
				rc = gzread(kcount_fp[i],kcounts+i,sizeof(uint32_t));
		}
		// choose first kmer
		size_t mindex = find_min(kmers,kcounts,batches);
		uint64_t distinct [nwords];
		memcpy(distinct,kmers + mindex*nwords, kmer_size);
		uint32_t count = kcounts[mindex];
		while(todo>0) {
			// read next kmer from kmer_fp[mindex]
			int rc = gzread(kmer_fp[mindex],kmers+mindex*nwords,kmer_size);
			if(rc == 0) {
				todo--;
				kcounts[mindex]=0;
				if (todo==0) break; // we're done
			}
			else if (rc == -1) {
				fprintf(stderr,"some error in gzread. mindex=%zi\n",mindex);
				todo--;
				kcounts[mindex]=0;
			}
			else
				gzread(kcount_fp[mindex],kcounts+mindex,sizeof(uint32_t));
			
			mindex = find_min(kmers,kcounts,batches);
			if (compare_kmers(kmers+mindex*nwords,distinct) == 0) // same kmer, update count
				count+= kcounts[mindex];
			else { // different kmer, output previous, setup next
				gzwrite(kmer_ofp,distinct,kmer_size);
				gzwrite(kcount_ofp,&count,sizeof(uint32_t));
				memcpy(distinct,kmers + mindex*nwords, kmer_size);
				count = kcounts[mindex];
			}
		}
		for(size_t i=0;i<batches; i++) {
			gzclose(kmer_fp[i]);
			gzclose(kcount_fp[i]);
		}
		gzclose(kmer_ofp);
		gzclose(kcount_ofp);
	}
}


inline uint64_t* canonicalize(uint64_t *packed, uint64_t *rcpack) {
	int cmp=0;
	for(size_t i=0;i<nwords;i++) {
		rcpack[i] = revcomp(packed[nwords-1-i]);
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

void serialize_batch(vector<uint64_t*> &kmer_buf,vector<size_t> &bin_tally,size_t batches,size_t k) {
	boost::thread_group tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j = (i+thread_bins>NBINS) ? NBINS : i+thread_bins;
		tg.create_thread(boost::bind(sortuniq,i,j,kmer_buf,bin_tally,batches,k));
	}
	tg.join_all();
	for(size_t i=0;i<NBINS;i++)
		bin_tally[i]=0;
}

void count_frequency(const size_t from, const size_t to, vector<uint32_t> *histogram, const size_t k) {
	for (size_t bin=from;bin<to;bin++) {
		// open the input file of counts
		// sort and uniqify the counts into the histogram vector of alternating counts and counts of counts
		gzFile kcount_fp;
		char kcount_file [100];
		sprintf(kcount_file,"%s/distinct_%zi-mer_counts.%zi.0",outprefix,k,bin);
		kcount_fp = gzopen(kcount_file, "rb");
		gzbuffer(kcount_fp,65536);
		uint32_t kcount;
		vector<uint32_t> counts;
		while (gzread(kcount_fp,&kcount,sizeof(uint32_t)) > 0)
			counts.push_back(kcount);
		gzclose(kcount_fp);
		sort(counts.begin(),counts.end());
		histogram[bin].push_back(counts[0]);
		histogram[bin].push_back(1);
		size_t j=0;
		for(size_t i=1;i<counts.size();i++) {
			if (counts[i] == histogram[bin][j])
				histogram[bin][j+1]++;
			else {
				histogram[bin].push_back(counts[i]);
				histogram[bin].push_back(1);
				j+=2;
			}
		}
	}
}

void merge_histograms(vector<uint32_t> *histogram) {
	// these histograms are very similar, so brute force is likely to be very good
	uint32_t k=1; // min value in any histogram
	uint32_t todo = NBINS;
	size_t offset [NBINS];
	size_t done [NBINS];
	for(size_t i=0;i<NBINS;i++) {
		offset[i]=0;
		done[i]=0;
	}
	while (todo>0) { // decrement whenever we reach end end of a histogram
		uint32_t v=0;
		for(size_t i=0;i<NBINS;i++) {
			if (done[i]==0) {
				// fprintf(stderr,"todo: %u i: %zi, done: %zi offset: %zi k: %u v: %u\n",todo,i,done[i],offset[i],k,v);
				if (histogram[i][offset[i]] == k) {
					v += histogram[i][offset[i]+1];
					offset[i] += 2;
					if (offset[i] == histogram[i].size()) {
						done[i]=1;
						todo--;
					}
				}
			}
		}
		// fprintf(stderr,"%u\t%u\t%u\n",k,v,todo);
		// if (k==3) {
		// 	exit(1);
		// }
		if (v>0) {
			printf("%u\t%u\n",k,v);
		}
		k++;
	}
}

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	// parse args
	if (argc != 6) {
		fprintf(stderr, "Usage: %s <input file> <k> <threads> <cap_bits> <output dir>\n", argv[0]);
		return 1;
	}
	gzFile fp;
	fp = gzopen(argv[1], "r");
	size_t k = atoi(argv[2]);
	threads = atoi(argv[3]);
	size_t cap_bits = atoi(argv[4]);
	outprefix = argv[5];
	// create output directory if it doesn't exist
	mkdir(outprefix,0755);

	thread_bins = NBINS/threads;
	
	nwords = ((k-1)>>5)+1; // number of packed 64bit words required to hold a k mer
	kmer_size = nwords*sizeof(uint64_t);
	size_t memory_cap = 1ULL<<cap_bits;
	size_t max_kmers_per_bin = memory_cap/kmer_size/NBINS;
	fprintf(stderr,"k: %zi, nwords: %zi, memory_cap: %zi, max_kmers_per_bin: %zi\n",k,nwords,memory_cap,max_kmers_per_bin);
	vector<size_t> bin_tally; // when a bin reaches max_kmers_per_bin, sort,uniq,count,deflate
	vector<uint64_t*> kmer_buf; // we'll use qsort from stdlib to sort kmers in place
	// reserve space
	kmer_buf.resize(NBINS);
	bin_tally.resize(NBINS);
	for(int i=0;i<NBINS;i++) {
		kmer_buf[i] = (uint64_t*) calloc(max_kmers_per_bin, kmer_size);
		bin_tally[i] = 0;
	}

	size_t batches=0;
	// process each seq from input
	kseq_t *seq = kseq_init(fp);
	int length;
	while ((length = kseq_read(seq)) >= 0) {
		if ((size_t)length >= k) {
			uint64_t packed [nwords];
			uint64_t rcpack [nwords];
			memset(packed,0,kmer_size);

			packed[0] = twobit(seq->seq.s[0]);
			for(size_t i=1; i<k; i++) { 
				register size_t word = i>>5;
				packed[word] <<= 2; // shift 2 bits to the left
				packed[word] |= twobit(seq->seq.s[i]); // insert the next two bits
			}

			// canonicalize the first kmer
			uint64_t *kmer = canonicalize(packed,rcpack);
			size_t bin = hashkmer(kmer,0); // kmer[0] & 255;

			memcpy(kmer_buf[bin] + nwords*bin_tally[bin],kmer,kmer_size);
			bin_tally[bin]++;
			// one more kmer would cause a seg fault, so check if we've filled up our allocated memory
			if (bin_tally[bin] == max_kmers_per_bin) {
				fprintf(stderr,"buffer is full in bin %zi\n",bin);
				batches++;
				serialize_batch(kmer_buf,bin_tally,batches,k);
				fprintf(stderr,"back from serialize_batch()\n");
			}

			// pack the rest of the sequence
			for(size_t i=k; i<(size_t)length;i++) {
				// shift each kmer by 2 bits
				packed[0] <<= 2;
				for(size_t word=1;word<nwords;word++) {
					packed[word-1] |= packed[word] >> 62;
					packed[word] <<= 2;
				}
				packed[nwords-1] |= twobit(seq->seq.s[i]);

				// canonicalize the ith kmer
				kmer = canonicalize(packed,rcpack);
				bin = hashkmer(kmer,0); // kmer[0] & 255;
				memcpy(kmer_buf[bin] + nwords*bin_tally[bin],kmer,kmer_size);
				bin_tally[bin]++;

				// check if we've filled up our allocated memory
				if (bin_tally[bin] == max_kmers_per_bin) {
					fprintf(stderr,"buffer is full in bin %zi\n",bin);
					batches++;
					serialize_batch(kmer_buf,bin_tally,batches,k);
					fprintf(stderr,"back from serialize_batch()\n");
				}
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	fprintf(stderr,"finished reading from input file\n");
	// process the last batch
	if (batches == 0) {
		serialize_batch(kmer_buf,bin_tally,batches,k);
	}
	else {
		batches++;
		serialize_batch(kmer_buf,bin_tally,batches,k);
		// move the merge out of here and into another program
		fprintf(stderr,"merging %zi batches\n",batches);
		// mergecounts(0,NBINS,batches,k)
		boost::thread_group merge_tg;
		for(size_t i=0;i<NBINS;i+=thread_bins) {
			size_t j = (i+thread_bins > NBINS) ? NBINS : i+thread_bins;
			merge_tg.create_thread(boost::bind(mergecounts,i,j,batches,k));
		}
		merge_tg.join_all();
	}
	fprintf(stderr,"starting stats\n");
	vector<uint32_t> histogram [NBINS];
	boost::thread_group stats_tg;
	for(size_t i=0;i<NBINS;i+=thread_bins) {
		size_t j = (i+thread_bins > NBINS) ? NBINS : i+thread_bins;
		stats_tg.create_thread(boost::bind(count_frequency,i,j,histogram,k));
	}
	stats_tg.join_all();
	// merge the histograms and output
	fprintf(stderr,"merge histogrmas\n");
	merge_histograms(histogram);
	return 0;
}