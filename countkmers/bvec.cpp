#include "bvec.h"

// constructor - given a sorted vector of distinct 64bit integers
bvec::bvec(vector<uint64_t>& vals) {
	uint64_t word_end = LITERAL_SIZE - 1;
	uint32_t word=0;
	uint64_t gap_words = vals.front()/LITERAL_SIZE;
	while (gap_words > FILLMASK) {
		words.push_back(BIT32 | FILLMASK);
		gap_words -= FILLMASK;
	}
	if (gap_words > 0) {
		words.push_back(gap_words | BIT32);
		word_end += gap_words*LITERAL_SIZE;
	} 
	for(vector<uint64_t>::iterator ii = vals.begin(); ii != vals.end(); ++ii) {
		if (*ii <= word_end) {
			word |= 1 << (uint32_t)(word_end - *ii);
		}
		else {
			if (word == ALL1S) {
				if ((words.back() & ONEFILL) && (words.back() != ONEFILL2 ))
					words.back()++;
				else
					words.push_back(ONEFILL1);
			}
			gap_words = (*ii - word_end)/LITERAL_SIZE;
			while (gap_words > FILLMASK) {
				words.push_back(BIT32 | FILLMASK);
				gap_words -= FILLMASK;
			}
			if (gap_words > 0)
				words.push_back((uint32_t)gap_words | BIT32);
			word_end += (gap_words+1)*LITERAL_SIZE;
			word = 1 << (uint32_t)(word_end - *ii);
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

// union of two bitvectors
bvec* bvec::union(bvec & bv) {
	bvec* res = new bvec;
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();
	// measure the first pair of words
	uint64_t a_pos = (*a & BIT32) ? *a & FILLMASK : 1;
	uint64_t b_pos = (*b & BIT32) ? *b & FILLMASK : 1;
	uint64_t res_pos=0;
	uint32_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	while(a != words.end() && b != bv.words.end()) {
		if (a_pos == b_pos) {
			if (*a & ONEFILL || *b & ONEFILL) {
				next_word = ONEFILL | (a_pos - res_pos);
			}
			else if (*a & BIT32) { // zero fill
				next_word = BIT32 | (a_pos - res_pos);
			}
			else { // literal word
				uint32_t u = *a | *b;
				next_word = (u == ALL1s) ? ONEFILL1 : u;
			}
			incr_a = true;
			incr_b = true;
			res_pos = a_pos;
		}
		else if (a_pos < b_pos) {
			if (*b & ONEFILL) {
				next_word = ONEFILL | (b_pos - res_pos);
				res_pos = b_pos;
				incr_a = true;
				incr_b = true;
			}
			else { // b is a 0-fill because it can't be a literal word and have b_pos > a_pos
				if (*a & ONEFILL) {
					next_word = ONEFILL | (a_pos - res_pos);
				}
				else if (*a & BIT32) { // a is 0-fill
					next_word = BIT32 | (a_pos - res_pos);
				}
				else { // literal
					next_word = *a;
				}
				res_pos = a_pos;
				incr_a = true;
			}
		}
		else { // a_pos > b_pos
			if (*a & ONEFILL) {
				next_word = ONEFILL | (a_pos - res_pos);
				res_pos = a_pos;
				incr_a = true;
				incr_b = true;
			}
			else { // a is a 0-fill because it can't be a literal word and have a_pos > b_pos
				if (*b & ONEFILL) {
					next_word = ONEFILL | (b_pos - res_pos);
				}
				else if (*b & BIT32) { // b is 0-fill
					next_word = BIT32 | (b_pos - res_pos);
				}
				else { // literal
					next_word = *b;
				}
				res_pos = b_pos;
				incr_b = true;
			}
		}
		if (incr_a) {
			while(a_pos <= res_pos && a != words.end()) {
				++a;
				if (a != words.end())
					a_pos += (*a & BIT32) ? *a & FILLMASK : 1;
			}
			incr_a = false;
		}
		if (incr_b) {
			while(b_pos <= res_pos && b != bv.words.end()) {
				++b;
				if (b != bv.words.end())
					b_pos += (*b & BIT32) ? *b & FILLMASK : 1;
			}
			incr_b = false;
		}
		if (next_word & BIT32
			&& res->words.size() > 0
			&& res->words.back() & BIT32
			&& res->words.back() & BIT31 == next_word & BIT31
			&& res->words.back() & FILLMASK < FILLMASK
		) {
			// merge fill words
			// but just don't exceed the capacity
			uint64_t n_words = res->words.back() & FILLMASK + next_word & FILLMASK;
			if (n_words > FILLMASK) {
				res->words.back() |= FILLMASK;
				n_words -= FILLMASK;
				next_word = (next_word & ONEFILL) | n_words;
			}
		}
		res->words.push_back(next_word);
	}
	if (a != words.end()) {
		if(*a & BIT32) {
			res->words.push_back((*a & ONEFILL) | (a_pos - res_pos));
			++a;
			if (a != words.end())
				res->words.insert(res->words.end(),a,words.end());
		}
		else
			res->words.insert(res->words.end(),a,words.end());
	}
	if (b != bv.words.end()) {
		if(*b & BIT32) {
			res->words.push_back((*b & ONEFILL) | (b_pos - res_pos));
			++b;
			if (b != bv.words.end())
				res->words.insert(res->words.end(),b,bv.words.end());
		}
		else
			res->words.insert(res->words.end(),b,bv.words.end());
	}
	return res;
}

bvec* bvec::intersect(bvec & bv) {
	bvec* res = new bvec;
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();
	// measure the first pair of words
	uint64_t a_pos = (*a & BIT32) ? *a & FILLMASK : 1;
	uint64_t b_pos = (*b & BIT32) ? *b & FILLMASK : 1;
	uint64_t res_pos=0;
	uint32_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	while(a != words.end() && b != bv.words.end()) {
		if (a_pos == b_pos) {
			if (*a & ONEFILL && *b & ONEFILL) {
				next_word = ONEFILL | (a_pos - res_pos);
			}
			else if (*a & BIT32) { // zero fill
				next_word = BIT32 | (a_pos - res_pos);
			}
			else { // literal word
				uint32_t u = *a & *b;
				next_word = (u == 0) ? BIT32 | 1 : u;
			}
			incr_a = true;
			incr_b = true;
			res_pos = a_pos;
		}
		else if (a_pos < b_pos) {
			if (*b & ONEFILL) {
				if(*a & ONEFILL) {
					next_word = ONEFILL | (a_pos - res_pos);
				}
				else if (*a & BIT32 == 0) {
					next_word = *a;
				}
				else {
					next_word = BIT32 | (a_pos - res_pos);
				}
				res_pos = a_pos;
				incr_b = true;
			}
			else { // b is a 0-fill because it can't be a literal word and have b_pos > a_pos
				next_word = BIT32 | (b_pos - res_pos);
				res_pos = b_pos;
				incr_a = true;
				incr_b = true;
			}
		}
		else { // a_pos > b_pos
			if (*a & ONEFILL) {
				if(*b & ONEFILL) {
					next_word = ONEFILL | (b_pos - res_pos);
				}
				else if (*b & BIT32 == 0) {
					next_word = *b;
				}
				else {
					next_word = BIT32 | (b_pos - res_pos);
				}
				res_pos = b_pos;
				incr_a = true;
			}
			else { // a is a 0-fill because it can't be a literal word and have a_pos > b_pos
				next_word = BIT32 | (a_pos - res_pos);
				res_pos = a_pos;
				incr_a = true;
				incr_b = true;
			}
		}
		if (incr_a) {
			while(a_pos <= res_pos && a != words.end()) {
				++a;
				if (a != words.end())
					a_pos += (*a & BIT32) ? *a & FILLMASK : 1;
			}
			incr_a = false;
		}
		if (incr_b) {
			while(b_pos <= res_pos && b != bv.words.end()) {
				++b;
				if (b != bv.words.end())
					b_pos += (*b & BIT32) ? *b & FILLMASK : 1;
			}
			incr_b = false;
			// ++b;
			// if (b != bv.words.end())
			// 	b_pos += (*b & BIT32) ? *b & FILLMASK : 1;
			// incr_b = false;
		}
		if (next_word & BIT32
			&& res->words.size() > 0
			&& res->words.back() & BIT32
			&& res->words.back() & BIT31 == next_word & BIT31
			&& res->words.back() & FILLMASK < FILLMASK
		) {
			// merge fill words
			// but just don't exceed the capacity
			uint64_t n_words = res->words.back() & FILLMASK + next_word & FILLMASK;
			if (n_words > FILLMASK) {
				res->words.back() |= FILLMASK;
				n_words -= FILLMASK;
				next_word = (next_word & ONEFILL) | n_words;
			}
		}
		res->words.push_back(next_word);
	}
	if (a != words.end()) {
		if(*a & BIT32) {
			res->words.push_back((*a & ONEFILL) | (a_pos - res_pos));
			++a;
			if (a != words.end())
				res->words.insert(res->words.end(),a,words.end());
		}
		else
			res->words.insert(res->words.end(),a,words.end());
	}
	if (b != bv.words.end()) {
		if(*b & BIT32) {
			res->words.push_back((*b & ONEFILL) | (b_pos - res_pos));
			++b;
			if (b != bv.words.end())
				res->words.insert(res->words.end(),b,bv.words.end());
		}
		else
			res->words.insert(res->words.end(),b,bv.words.end());
	}
	return res;
}