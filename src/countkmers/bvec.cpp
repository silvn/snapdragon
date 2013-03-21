#include "bvec.h"

// constructor - given a sorted vector of distinct 32bit integers
bvec::bvec(vector<uint32_t>& vals) {
	uint32_t word_end = LITERAL_SIZE - 1;
	uint32_t word=0;
	uint32_t gap_words = vals.front()/LITERAL_SIZE;
	// while (gap_words > FILLMASK) {
	// 	words.push_back(BIT32 | FILLMASK);
	// 	gap_words -= FILLMASK;
	// }
	if (gap_words > 0) {
		words.push_back(gap_words | BIT32);
		word_end += gap_words*LITERAL_SIZE;
	} 
	for(vector<uint32_t>::iterator ii = vals.begin(); ii != vals.end(); ++ii) {
		if (*ii <= word_end)
			word |= 1 << (word_end - *ii);
		else {
			if (word == ALL1S)
				if ((words.back() & ONEFILL) && (words.back() != ONEFILL2 ))
					words.back()++;
				else
					words.push_back(ONEFILL1);
			else
				words.push_back(word);
			gap_words = (*ii - word_end)/LITERAL_SIZE;
			if (gap_words > 0)
				words.push_back(gap_words | BIT32);
			word_end += (gap_words+1)*LITERAL_SIZE;
			word = 1 << (word_end - *ii);
		}
	}
	// add the last word
	if (word == ALL1S) {
		if ((words.size() != 0) && ((words.back() & ONEFILL) == ONEFILL) && (words.back() != ONEFILL2 ))
			words.back()++;
		else
			words.push_back(ONEFILL1);
	}
	else {
		words.push_back(word);
	}
	count = vals.size();
	size = word_end+1;
}
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
		if (*ii <= word_end)
			word |= 1 << (uint32_t)(word_end - *ii);
		else {
			if (word == ALL1S)
				if ((words.back() & ONEFILL) && (words.back() != ONEFILL2 ))
					words.back()++;
				else
					words.push_back(ONEFILL1);
			else
				words.push_back(word);
			gap_words = (*ii - word_end)/LITERAL_SIZE;
			while (gap_words > FILLMASK) {
				words.push_back(BIT32 | FILLMASK);
				gap_words -= FILLMASK;
			}
			if (gap_words > 0)
				words.push_back(BIT32 | (uint32_t)gap_words);
			word_end += (gap_words+1)*LITERAL_SIZE;
			word = 1 << (uint32_t)(word_end - *ii);
		}
	}
	// add the last word
	if (word == ALL1S) {
		if ((words.size() != 0) && ((words.back() & ONEFILL) == ONEFILL) && (words.back() != ONEFILL2 ))
			words.back()++;
		else
			words.push_back(ONEFILL1);
	}
	else {
		words.push_back(word);
	}
	count = vals.size();
	size = word_end+1;
}

void bvec::matchSize(bvec &bv) {
	if (size < bv.size) {
		uint64_t gap_words = (bv.size - size)/LITERAL_SIZE;
		while (gap_words > FILLMASK) {
			words.push_back(BIT32 | FILLMASK);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			words.push_back(BIT32 | (uint32_t)gap_words);
		size = bv.size;
	}
	else if (size > bv.size) {
		uint64_t gap_words = (size - bv.size)/LITERAL_SIZE;
		while (gap_words > FILLMASK) {
			bv.words.push_back(BIT32 | FILLMASK);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			bv.words.push_back(BIT32 | (uint32_t)gap_words);
		bv.size = size;
	}
}

// in place version of the bitwise OR operator.
void bvec::operator|=(bvec& bv) {
	// ensure that both bvecs are the same size
	this->matchSize(bv);
	if (size == 0)
		return;
	
	vector<uint32_t> res; // fill this then swap with this.words
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();

	// maintain the end position of the current word
	uint64_t a_pos = (*a & BIT32) ? *a & FILLMASK : 1;
	uint64_t b_pos = (*b & BIT32) ? *b & FILLMASK : 1;
	uint64_t res_pos=0;
	uint32_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	uint64_t last_pos = size/LITERAL_SIZE;
	while(res_pos != last_pos) {
		if (incr_a) {
			while(a_pos <= res_pos) {
				++a;
				a_pos += (*a & BIT32) ? *a & FILLMASK : 1;
			}
			incr_a = false;
		}
		if (incr_b) {
			while(b_pos <= res_pos) {
				++b;
				b_pos += (*b & BIT32) ? *b & FILLMASK : 1;
			}
			incr_b = false;
		}
		if (a_pos == b_pos) {
			if ((*a & ONEFILL) == ONEFILL || (*b & ONEFILL) == ONEFILL)
				next_word = ONEFILL | (a_pos - res_pos);
			else
				if (*a & BIT32)
					if (*b & BIT32) // zero fill
						next_word = BIT32 | (a_pos - res_pos);
					else
						next_word = *b;
				else
					if (*b & BIT32) // zero fill
						next_word = *a;
					else {
						uint32_t u = *a | *b;
						next_word = (u == ALL1S) ? ONEFILL1 : u;
					}
			incr_a = true;
			incr_b = true;
			res_pos = a_pos;
		}
		else if (a_pos < b_pos) {
			if ((*b & ONEFILL) == ONEFILL) {
				next_word = ONEFILL | (b_pos - res_pos);
				res_pos = b_pos;
				incr_a = true;
				incr_b = true;
			}
			else { // b is a 0-fill or a literal word
				if ((*a & ONEFILL) == ONEFILL)
					next_word = ONEFILL | (a_pos - res_pos);
				else if (*a & BIT32) // a is 0-fill
					next_word = BIT32 | (a_pos - res_pos);
				else // literal
					next_word = *a;
				res_pos = a_pos;
				incr_a = true;
			}
		}
		else { // a_pos > b_pos
			if ((*a & ONEFILL) == ONEFILL) {
				next_word = ONEFILL | (a_pos - res_pos);
				res_pos = a_pos;
				incr_a = true;
				incr_b = true;
			}
			else { // a is a 0-fill or a literal word
				if ((*b & ONEFILL) == ONEFILL)
					next_word = ONEFILL | (b_pos - res_pos);
				else if (*b & BIT32) // b is 0-fill
					next_word = BIT32 | (b_pos - res_pos);
				else // literal
					next_word = *b;
				res_pos = b_pos;
				incr_b = true;
			}
		}
		if ((next_word & BIT32)
			&& res.size() > 0
			&& (res.back() & BIT32)
			&& (res.back() & BIT31) == (next_word & BIT31)
			&& (res.back() & FILLMASK) < FILLMASK
		) {
			// merge fill words
			// but just don't exceed the capacity
			uint64_t n_words = (res.back() & FILLMASK) + (next_word & FILLMASK);
			if (n_words >= FILLMASK) {
				res.back() |= FILLMASK;
				n_words -= FILLMASK;
				if (n_words > 0)
					 res.push_back((next_word & ONEFILL) | n_words);
			}
			else
				res.push_back(next_word);
		}
		else
			res.push_back(next_word);
	}
	words.swap(res);
	count=0;
}

// in place version of the bitwise AND operator.
void bvec::operator&=(bvec& bv) {
	// ensure that both bvecs are the same size
	this->matchSize(bv);
	if (size == 0)
		return;
	
	vector<uint32_t> res; // fill this then swap with this.words
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();

	// maintain the end position of the current word
	uint64_t a_pos = (*a & BIT32) ? *a & FILLMASK : 1;
	uint64_t b_pos = (*b & BIT32) ? *b & FILLMASK : 1;
	uint64_t res_pos=0;
	uint32_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	uint64_t last_pos = size/LITERAL_SIZE;
	while(res_pos != last_pos) {
		if (incr_a) {
			while(a_pos <= res_pos) {
				++a;
				a_pos += (*a & BIT32) ? *a & FILLMASK : 1;
			}
			incr_a = false;
		}
		if (incr_b) {
			while(b_pos <= res_pos) {
				++b;
				b_pos += (*b & BIT32) ? *b & FILLMASK : 1;
			}
			incr_b = false;
		}
		if (a_pos == b_pos) {
			if ((*a & ONEFILL) == ONEFILL && (*b & ONEFILL) == ONEFILL)
				next_word = ONEFILL | (a_pos - res_pos);
			else if ((*a & BIT32) || (*b & BIT32))
				next_word = BIT32 | (a_pos - res_pos);
			else {
				uint32_t u = *a & *b;
				next_word = (u == 0) ? BIT32 | 1 : u;
			}
			incr_a = true;
			incr_b = true;
			res_pos = a_pos;
		}
		else if (a_pos < b_pos) {
			if ((*b & ONEFILL) == ONEFILL) {
				if((*a & ONEFILL) == ONEFILL)
					next_word = ONEFILL | (a_pos - res_pos);
				else if (*a & BIT32)
					next_word = BIT32 | (a_pos - res_pos);
				else
					next_word = *a;
				res_pos = a_pos;
				incr_a = true;
			}
			else { // b is a 0-fill because it can't be a literal word and have b_pos > a_pos
				next_word = BIT32 | (b_pos - res_pos);
				res_pos = b_pos;
				incr_a = true;
				incr_b = true;
			}
		}
		else { // a_pos > b_pos
			if ((*a & ONEFILL) == ONEFILL) {
				if((*b & ONEFILL) == ONEFILL)
					next_word = ONEFILL | (b_pos - res_pos);
				else if (*b & BIT32)
					next_word = BIT32 | (b_pos - res_pos);
				else
					next_word = *b;
				res_pos = b_pos;
				incr_b = true;
			}
			else { // a is a 0-fill because it can't be a literal word and have a_pos > b_pos
				next_word = BIT32 | (a_pos - res_pos);
				res_pos = a_pos;
				incr_a = true;
				incr_b = true;
			}
		}
		if ((next_word & BIT32)
			&& res.size() > 0
			&& (res.back() & BIT32)
			&& (res.back() & BIT31) == (next_word & BIT31)
			&& (res.back() & FILLMASK) < FILLMASK
		) {
			// merge fill words
			// but just don't exceed the capacity
			uint64_t n_words = (res.back() & FILLMASK) + (next_word & FILLMASK);
			if (n_words >= FILLMASK) {
				res.back() |= FILLMASK;
				n_words -= FILLMASK;
				if (n_words > 0)
					 res.push_back((next_word & ONEFILL) | n_words);
			}
			else
				res.push_back(next_word);
		}
		else
			res.push_back(next_word);
	}
	words.swap(res);
	count=0;
}

bvec* bvec::operator|(bvec& rhs) {
	bvec *res = new bvec();
	if (words.size() > rhs.words.size()) {
		res->copy(rhs);
		*res |= *this;
		return res;
	}
	res->copy(*this);
	*res |= rhs;
	return res;
}
bvec* bvec::operator&(bvec& rhs) {
	bvec *res = new bvec();
	if (words.size() > rhs.words.size()) {
		res->copy(rhs);
		*res &= *this;
		return res;
	}
	res->copy(*this);
	*res &= rhs;
	return res;
}