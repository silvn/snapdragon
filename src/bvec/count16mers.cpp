#define NBINS 256
#define NBITS 56
#define NTHREADS 16
#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <boost/thread.hpp>
#include "kseq.h"
using namespace std;

vector<uint32_t> kmers [NBINS]; // use the first 8 bits to choose a bin

vector<uint32_t> n_unique;

void processBins(int from, int to) {
	for(int i=from; i<to; i++) {
		sort(kmers[i].begin(),kmers[i].end());
		if (kmers[i].size() > 0) {
			uint32_t tally=1;
			uint32_t mer = kmers[i][0];
			for(int j=1; j<kmers[i].size(); j++) {
				if (mer == kmers[i][j]) {
					tally++;
				}
				else {
					if (tally == 1) n_unique[i]++;
//					printf("%llu %u\n",mer,tally);
					mer = kmers[i][j];
					tally=1;
				}
			}
		}
	}
}

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	vector<uint32_t> twoBit;
	twoBit.resize(256,0);
	twoBit[99] = 1;
	twoBit[103] = 2;
	twoBit[116] = 3;
	twoBit[67] = 1;
	twoBit[71] = 2;
	twoBit[84] = 3;
	gzFile fp;
	kseq_t *seq;
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		if (seq->seq.l>=16) {
			uint64_t mer = twoBit[seq->seq.s[0]];
			for(int i=1;i<16;i++) { // fill first 32bit word
				mer <<= 2; // shift 2 bits to the left
				mer |= twoBit[seq->seq.s[i]]; // insert the next two bits
			}
			size_t bin = mer >> NBITS;
			kmers[bin].push_back(mer);
			for(int i=16;i<seq->seq.l;i++) {
				mer <<= 2;
				mer |= twoBit[seq->seq.s[i]];
				bin = mer >> NBITS;
				kmers[bin].push_back(mer);
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	// sort each bin and count 16mers
	n_unique.resize(NBINS,0);
	boost::thread_group tg;
	int groupsize = NBINS/NTHREADS;
	for(int i=0; i<NBINS; i+=groupsize) {
		int j = i+groupsize;
		tg.create_thread(boost::bind(processBins, i, j>NBINS ? j=NBINS : j));
	}
	tg.join_all();
	uint32_t total=0;
	for(int i=0; i<NBINS; i++) {
		total += n_unique[i];
	}
	printf("number of unique 16mers: %i\n",total);
	return 0;
}
