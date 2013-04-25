#define NBINS 256
#define NBITS 56
#include <zlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include "kseq.h"
#include "bvec32.h"
using namespace std;

inline uint32_t revcomp(uint32_t val) {
	// bitwise reverse complement of values from 0 to 255
	static const uint32_t table[256] = {
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
		240,112,176,48,208,80,144,16,224,96,160,32,192,64,128,0};
	return (table[val&0xFFUL]<<24) | (table[(val>>8)&0xFFUL]<<16) |
		(table[(val>>16)&0xFFUL]<<8) | (table[(val>>24)&0xFFUL]);
}

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
	unsigned int k;
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <in.seq> <mer> <qual>\n", argv[0]);
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
			uint32_t mer = twoBit[seq->seq.s[0]];
			for(int i=1; i<k; i++) { // fill first word
				mer <<= 2; // shift 2 bits to the left
				mer |= twoBit[seq->seq.s[i]]; // insert the next two bits
			}
			vector<uint32_t> kmers;
			kmers.reserve(length - k + 1);
			uint32_t rem = revcomp(mer);
			kmers.push_back(mer < rem ? mer : rem);
			for(int i=k; i<length; i++) {
				mer <<= 2;
				mer |= twoBit[seq->seq.s[i]];
				mer &= kmask;
				rem = revcomp(mer);
				kmers.push_back(mer < rem ? mer : rem);
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

	// 
	// vector<bvec32*>::iterator si = seq2vec.begin();
	// bvec32 *merged = *si;
	// ++si;
	// while (si != seq2vec.end()) {
	// 	*merged |= **si;
	// 	++si;
	// }
	// printf("finished ORing all the bvecs. %u distinct %i-mers\n",merged->cnt(),k);
	return 0;
}
