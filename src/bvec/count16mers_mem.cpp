#define NBINS 256
#define NBITS 24
#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/thread.hpp>
#include "kseq.h"
using namespace std;

vector<uint32_t> kmers [NBINS]; // use the first 8 bits to choose a bin
vector<uint32_t> kmer_keys [NBINS];
vector<uint32_t> kmer_counts [NBINS];

void processBins(int from, int to) {
	for(int i=from; i<to; i++) {
		sort(kmers[i].begin(),kmers[i].end());
		vector<uint32_t> keys;
		vector<uint32_t> counts;
		if (kmers[i].size() > 0) {
			uint32_t tally=1;
			uint32_t mer = kmers[i][0];
			for(int j=1; j<kmers[i].size(); j++) {
				if (mer == kmers[i][j]) {
					tally++;
				}
				else {
					keys.push_back(mer);
					counts.push_back(tally);
					mer = kmers[i][j];
					tally=1;
				}
			}
			kmers[i].clear();
			vector<uint32_t>::iterator ki = keys.begin();
			vector<uint32_t>::iterator ci = counts.begin();
			vector<uint32_t>::iterator kki = kmer_keys[i].begin();
			vector<uint32_t>::iterator kci = kmer_counts[i].begin();
			vector<uint32_t> merged_keys;
			vector<uint32_t> merged_counts;
			while(ki != keys.end() && kki != kmer_keys[i].end()) {
				if (*ki < *kki) {
					merged_keys.push_back(*ki);
					merged_counts.push_back(*ci);
					*ki++;*ci++;
				}
				else if(*kki < *ki) {
					merged_keys.push_back(*kki);
					merged_counts.push_back(*kci);
					*kki++;*kci++;
				}
				else {
					merged_keys.push_back(*ki);
					merged_counts.push_back(*ci+*kci);
					*ki++;*ci++;*kki++;*kci++;
				}
			}
			if (ki != keys.end()) {
				merged_keys.insert(merged_keys.end(),ki,keys.end());
				merged_counts.insert(merged_counts.end(),ci,counts.end());
			}
			if (kki != kmer_keys[i].end()) {
				merged_keys.insert(merged_keys.end(),kki,kmer_keys[i].end());
				merged_counts.insert(merged_counts.end(),kci,kmer_counts[i].end());
			}
			merged_keys.swap(kmer_keys[i]);
			merged_counts.swap(kmer_counts[i]);
		}
	}
}

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	if (argc != 4) {
		fprintf(stderr, "Usage: %s <in.seq> <threads> <max kmers before uniqify>\n", argv[0]);
		return 1;
	}

	// setup 2bit vector
	vector<uint32_t> twoBit;
	twoBit.resize(256,0);
	twoBit[99] = 1;  // c
	twoBit[103] = 2; // g
	twoBit[116] = 3; // t
	twoBit[67] = 1;  // C
	twoBit[71] = 2;  // G
	twoBit[84] = 3;  // T


	gzFile fp;
	kseq_t *seq;
	int length;
	fp = gzopen(argv[1], "r");
	uint32_t nthreads = atoi(argv[2]);
	int groupsize = NBINS/nthreads;
	uint32_t max_unsorted = atoi(argv[3]);
	// reserve space in kmers vector
	for(int i=0;i<NBINS;i++) {
		kmers[i].reserve(max_unsorted);
	}

	seq = kseq_init(fp);
	while ((length = kseq_read(seq)) >= 0) {
		if (length >= 16) {
			uint32_t mer = twoBit[seq->seq.s[0]];
			for(int i=1; i<16; i++) { // fill first 32bit word
				mer <<= 2; // shift 2 bits to the left
				mer |= twoBit[seq->seq.s[i]]; // insert the next two bits
			}
			size_t bin = mer >> NBITS;
			kmers[bin].push_back(mer);
			for(int i=16; i<length; i++) {
				if (kmers[bin].size() == max_unsorted) {
					// once you reach this cutoff, process all bins
					// in multiple threads and merge the kmer tallies
					boost::thread_group tg;
					for(int i=0; i<NBINS; i+=groupsize) {
						int j = i+groupsize;
						tg.create_thread(boost::bind(processBins, i, j>NBINS ? j=NBINS : j));
					}
					tg.join_all();
				}
				mer <<= 2;
				mer |= twoBit[seq->seq.s[i]];
				bin = mer >> NBITS;
				kmers[bin].push_back(mer);
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	boost::thread_group tg;
	for(int i=0; i<NBINS; i+=groupsize) {
		int j = i+groupsize;
		tg.create_thread(boost::bind(processBins, i, j>NBINS ? j=NBINS : j));
	}
	tg.join_all();
	printf("Finished counting kmers\n");
	return 0;
}
