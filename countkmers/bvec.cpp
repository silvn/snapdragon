#include "bvec.h"

// constructor - given a sorted vector of distinct 64bit integers
bvec::bvec(vector<uint64_t>& vals) : count(0) {
	uint64_t word_end = LITERAL_SIZE - 1;
	uint32_t word=0;
	uint32_t gap_words = vals.front()/LITERAL_SIZE;
	if (gap_words > 0) {
		words.push_back(gap_words | BIT32);
		word_end += gap_words*LITERAL_SIZE;
	} 
	for(vector<uint64_t>::iterator ii = vals.begin(); ii != vals.end(); ++ii) {
		if (*ii <= word_end) {
			word |= 1 << (word_end - *ii);
		}
		else {
			if (word == ALL1S) {
				if (words.back() & ONEFILL)
					words.back()++;
				else
					words.push_back(ONEFILL1);
			}
			gap_words = (*ii - word_end)/LITERAL_SIZE;
			if (gap_words > 0)
				words.push_back(gap_words | BIT32);
			word_end += (gap_words+1)*LITERAL_SIZE;
			word = 1 << (word_end - *ii);
		}
	}
	count = vals.size();
}

// count the number of set bits
uint64_t bvec::cnt() {
	if (count == 0) 
		for(vector<uint32_t>::iterator it = words.begin(); it != words.end(); ++it)
			count += (*it & BIT32) ? (*it & BIT31) ? (*it & FILLMASK) * LITERAL_SIZE : 0 : __builtin_popcount(*it) - 1;
	return count;
}
