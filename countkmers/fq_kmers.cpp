#define NBINS 256
#define NBITS 56
#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include "kseq.h"
#include "bvec32.h"
using namespace std;


KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	unsigned int k;
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <in.seq> <mer> <qual>\n", argv[0]);
		return 1;
	}

	// setup 2bit vector
	vector<uint32_t> twobit;
	twobit.resize(256,0);
	twobit[99] = 1;  // c
	twobit[103] = 2; // g
	twobit[116] = 3; // t
	twobit[67] = 1;  // C
	twobit[71] = 2;  // G
	twobit[84] = 3;  // T

	gzFile fp;
	kseq_t *seq;
	int length;
	fp = gzopen(argv[1], "r");
	k = atoi(argv[2]);
	if (k > 16) {
		fprintf(stderr, "this program only works on up to 16mers\n");
		return 1;
	}
	unsigned int qual = atoi(argv[3]);
	seq = kseq_init(fp);
	vector<bvec32*> seq2vec;
	uint32_t kmask = (k==16) ? 0xFFFFFFFF : (1 << (2*k)) - 1;
	vector<uint32_t> allkmers;
	while ((length = kseq_read(seq)) >= 0) {
		if (length >= k) {
			uint32_t mer = twobit[seq->seq.s[0]];
			for(int i=1; i<k; i++) { // fill first word
				mer <<= 2; // shift 2 bits to the left
				mer |= twobit[seq->seq.s[i]]; // insert the next two bits
			}
			vector<uint32_t> kmers;
			kmers.reserve(length - k + 1);
			kmers.push_back(mer);
			for(int i=k; i<length; i++) {
				mer <<= 2;
				mer |= twobit[seq->seq.s[i]];
				mer &= kmask;
				kmers.push_back(mer);
			}
			sort(kmers.begin(),kmers.end());
			vector<uint32_t> merged;
			merged.reserve(length - k + 1);
			vector<uint32_t>::iterator ki = kmers.begin();
			merged.push_back(*ki);
			ki++;
			while (ki != kmers.end()) {
				if (*ki != merged.back())
					merged.push_back(*ki);
				ki++;
			}
			seq2vec.push_back(new bvec32(merged));
			allkmers.insert(allkmers.end(),merged.begin(),merged.end());
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	printf("finished reading %i-mers from %zi sequences\n",k,seq2vec.size());
	sort(allkmers.begin(),allkmers.end());
	printf("sorted %zi allkmers\n",allkmers.size());
	vector<uint32_t> allmerged;
	vector<uint32_t>::iterator aki = allkmers.begin();
	allmerged.push_back(*aki);
	while(aki != allkmers.end()) {
		if (*aki != allmerged.back())
			allmerged.push_back(*aki);
		aki++;
	}
	printf("after merging via sorting %zi distinct kmers\n",allmerged.size());
	vector<bvec32*>::iterator si = seq2vec.begin();
	bvec32 *merged = *si;
	++si;
	while (si != seq2vec.end()) {
		*merged |= **si;
		++si;
	}
	printf("finished ORing all the bvecs. %u distinct %i-mers\n",merged->cnt(),k);
	return 0;
}
