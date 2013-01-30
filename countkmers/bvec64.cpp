#include "bvec64.h"

// constructor - given a sorted vector of distinct 64bit integers
bvec::bvec(vector<uint64_t>& vals) {
	count = vals.size();
	// if the density is too low, run length encoding will take MORE space
	if (low_density(vals)) {
		if (DEBUG) printf("constructed non rle\n");
		words = vals;
		rle = false;
		size = 0; // irrelevant?
	}
	else {
		if (DEBUG) printf("construct_rle\n");
		construct_rle(vals);
	}
}

bool bvec::low_density(vector<uint64_t>& vals) {
	if (DEBUG) printf("low_density() %u/%u %c %f\n", (uint64_t)vals.size(),
		vals.back() - vals.front() + 1,
			((double)vals.size()/(double)(vals.back() - vals.front() + 1) < 1.0/(double)LITERAL_SIZE) ? '<' : '>',
				1.0/(double)LITERAL_SIZE);
	return (double)vals.size()/(double)(vals.back() - vals.front() + 1) < 1.0/(double)LITERAL_SIZE;
	
}

void bvec::print() {
	printf("rle: %c\n",rle ? 'T' : 'F');
	printf("words:\n");
	for(int i=0;i<words.size();i++) {
		printf(" %i ",i);
		if (rle)
			if ((words[i] & ONEFILL) == ONEFILL) 
				printf("1-fill %u %u\n",words[i] & FILLMASK, LITERAL_SIZE*(words[i] & FILLMASK));
			else if (words[i] & BIT1)
				printf("0-fill %u %u\n",words[i] & FILLMASK, LITERAL_SIZE*(words[i] & FILLMASK));
			else
				printf("literal %u\n",words[i]);
		else
			printf("direct %u\n",words[i]);
	}
}

void bvec::construct_rle(vector<uint64_t>& vals) {
	rle = true;
	uint64_t word_end = LITERAL_SIZE - 1;
	uint64_t word=0;
	uint64_t gap_words = vals.front()/LITERAL_SIZE;
	if (gap_words > 0) {
		word_end += LITERAL_SIZE*gap_words;
		while (gap_words > FILLMASK) {
			words.push_back(ZEROFULL);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			words.push_back(gap_words | BIT1);
	}
	for(vector<uint64_t>::iterator ii = vals.begin(); ii != vals.end(); ++ii) {
		if (*ii == word_end)
			word |= 1;
		else if (*ii < word_end)
			word |= ((uint64_t)1 << (word_end - *ii));
		else {
			if (word == ALL1S)
				if ((words.size() != 0) &&
                    ((words.back() & ONEFILL) == ONEFILL) &&
                    (words.back() != ONEFULL))
					words.back()++;
				else
					words.push_back(ONEFILL1);
			else
				words.push_back(word);
			gap_words = (*ii - word_end - 1)/LITERAL_SIZE;
			while (gap_words > FILLMASK) {
				words.push_back(ZEROFULL);
				gap_words -= FILLMASK;
			}
			if (gap_words > 0)
				words.push_back(gap_words | BIT1);
			word_end += (gap_words+1)*LITERAL_SIZE;
			word = (word_end - *ii == LITERAL_SIZE)
                ? 1 : (uint64_t)1 << (word_end - *ii);
		}
	}
	// add the last word
	if (word == ALL1S) {
		if ((words.size() != 0) &&
            ((words.back() & ONEFILL) == ONEFILL) &&
                (words.back() != ONEFULL))
			words.back()++;
		else
			words.push_back(ONEFILL1);
	} else {
		words.push_back(word);
	}
	size = word_end+1;
}

void bvec::compress() {
    if (rle) { /* Throw exception? */ return; }
	vector<uint64_t> tmp;
	tmp.swap(words);
	construct_rle(tmp);
    rle = true;
}

vector<uint64_t>& bvec::get_words() {
    return words;
}

void bvec::decompress() {
    if (!rle) { /* Throw exception? */ return; }
	// retrieve the set bits from the compressed vector
	vector<uint64_t> res;
	res.reserve(cnt());
	uint64_t pos=0;
	for(vector<uint64_t>::iterator ii = words.begin(); ii != words.end(); ++ii) {
		if ((*ii & ONEFILL) == ONEFILL) {
			uint64_t n_ones = LITERAL_SIZE*(*ii & FILLMASK);
			for(uint64_t i=0;i < n_ones; i++) {
				res.push_back(pos);
				pos++;
			}
		}
		else if (*ii & BIT1)
			pos += LITERAL_SIZE*(*ii & FILLMASK);
		else { // desconstruct the literal word
			for(uint64_t bit=1; bit<=LITERAL_SIZE; bit++)
				if (*ii & ((uint64_t)1 << (LITERAL_SIZE-bit)))
					res.push_back(pos+bit-1);
			pos += LITERAL_SIZE;
		}
	}
	words.swap(res);
	count = words.size();
	rle = false;
}

// in place version of the bitwise OR operator.
void bvec::operator|=(bvec& bv) {
	// decide which version we'll be using
	if (rle)
		if (bv.rle)
			rle_OR_rle(bv);
		else
			rle_OR_non(bv);
	else
		if (bv.rle)
			non_OR_rle(bv);
		else
			non_OR_non(bv);
}
// in place version of the bitwise AND operator.
void bvec::operator&=(bvec& bv) {
	// decide which version we'll be using
	if (rle)
		if (bv.rle)
			rle_AND_rle(bv);
		else
			rle_AND_non(bv);
	else
		if (bv.rle)
			non_AND_rle(bv);
		else
			non_AND_non(bv);
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

void bvec::non_OR_non(bvec& bv) {
	
	vector<uint64_t> res;
	vector<uint64_t>::iterator a = words.begin();
	vector<uint64_t>::iterator b = bv.words.begin();
	res.push_back(*a < *b ? *a : *b);
	while(a != words.end() && b != bv.words.end()) {
		if (*a < *b) {
			if (*a != res.back())
				res.push_back(*a);
			++a;
		}
		else if (*b < *a) {
			if (*b != res.back())
				res.push_back(*b);
			++b;
		}
		else {
			if (*a != res.back())
				res.push_back(*a);
			++a;
			++b;
		}
	}
	if (a != words.end())
		res.insert(res.end(),a,words.end());
	else if (b != bv.words.end())
		res.insert(res.end(),b,bv.words.end());

	count = res.size();
	// TODO: check if it's worth compressing

	words.swap(res);
}

void bvec::non_AND_non(bvec& bv) {
	
	vector<uint64_t> res;
	vector<uint64_t>::iterator a = words.begin();
	vector<uint64_t>::iterator b = bv.words.begin();
	res.push_back(*a < *b ? *a : *b);
	while(a != words.end() && b != bv.words.end()) {
		if (*a < *b)
			++a;
		else if (*b < *a)
			++b;
		else {
			res.push_back(*a);
			++a;
			++b;
		}
	}
	count = res.size();
	words.swap(res);
}

void bvec::non_AND_rle(bvec& bv) {
	// decompress
	// run non_AND_non
	bvec *tmp = new bvec();
	tmp->copy(bv);
	tmp->decompress();
	non_AND_non(*tmp);
	delete tmp;
}
void bvec::rle_AND_non(bvec& bv) {
	// decompress
	// run non_AND_non
	decompress();
	non_AND_non(bv);
}
void bvec::non_OR_rle(bvec& bv) {
	// compress
	// run rle_OR_rle
	compress();
	rle_OR_rle(bv);
}
void bvec::rle_OR_non(bvec& bv) {
	// compress
	// run rle_OR_rle
	bvec *tmp = new bvec();
	tmp->copy(bv);
	tmp->compress();
	rle_OR_rle(*tmp);
	delete tmp;
}

void bvec::matchSize(bvec &bv) {
	if (size < bv.size) {
		uint64_t gap_words = (bv.size - size)/LITERAL_SIZE;
		while (gap_words > FILLMASK) {
			words.push_back(ZEROFULL);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			words.push_back(BIT1 | gap_words);
		size = bv.size;
	}
	else if (size > bv.size) {
		uint64_t gap_words = (size - bv.size)/LITERAL_SIZE;
		while (gap_words > FILLMASK) {
			bv.words.push_back(ZEROFULL);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			bv.words.push_back(BIT1 | gap_words);
		bv.size = size;
	}
}

void bvec::rle_OR_rle(bvec& bv) {
	// ensure that both bvecs are the same size
	this->matchSize(bv);
	if (size == 0)
		return;
	
	vector<uint64_t> res; // fill this then swap with this.words
	vector<uint64_t>::iterator a = words.begin();
	vector<uint64_t>::iterator b = bv.words.begin();

	// maintain the end position of the current word
	uint64_t a_pos = (*a & BIT1) ? *a & FILLMASK : 1;
	uint64_t b_pos = (*b & BIT1) ? *b & FILLMASK : 1;
	uint64_t res_pos=0;
	uint64_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	uint64_t last_pos = size/LITERAL_SIZE;
	while(res_pos != last_pos) {
		if (incr_a) {
			while(a_pos <= res_pos) {
				++a;
				a_pos += (*a & BIT1) ? *a & FILLMASK : 1;
			}
			incr_a = false;
		}
		if (incr_b) {
			while(b_pos <= res_pos) {
				++b;
				b_pos += (*b & BIT1) ? *b & FILLMASK : 1;
			}
			incr_b = false;
		}
		if (a_pos == b_pos) {
			if ((*a & ONEFILL) == ONEFILL || (*b & ONEFILL) == ONEFILL)
				next_word = ONEFILL | (a_pos - res_pos);
			else
				if (*a & BIT1)
					if (*b & BIT1) // zero fill
						next_word = BIT1 | (a_pos - res_pos);
					else
						next_word = *b;
				else
					if (*b & BIT1) // zero fill
						next_word = *a;
					else {
						uint64_t u = *a | *b;
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
				else if (*a & BIT1) // a is 0-fill
					next_word = BIT1 | (a_pos - res_pos);
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
				else if (*b & BIT1) // b is 0-fill
					next_word = BIT1 | (b_pos - res_pos);
				else // literal
					next_word = *b;
				res_pos = b_pos;
				incr_b = true;
			}
		}
		if ((next_word & BIT1)
			&& res.size() > 0
			&& (res.back() & BIT1)
			&& (res.back() & BIT2) == (next_word & BIT2)
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
	
	// decide whether to decompress
}

// in place version of the bitwise AND operator.
void bvec::rle_AND_rle(bvec& bv) {
	// ensure that both bvecs are the same size
	this->matchSize(bv);
	if (size == 0)
		return;
	
	vector<uint64_t> res; // fill this then swap with this.words
	vector<uint64_t>::iterator a = words.begin();
	vector<uint64_t>::iterator b = bv.words.begin();

	// maintain the end position of the current word
	uint64_t a_pos = (*a & BIT1) ? *a & FILLMASK : 1;
	uint64_t b_pos = (*b & BIT1) ? *b & FILLMASK : 1;
	uint64_t res_pos=0;
	uint64_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	uint64_t last_pos = size/LITERAL_SIZE;
	while(res_pos != last_pos) {
		if (incr_a) {
			while(a_pos <= res_pos) {
				++a;
				a_pos += (*a & BIT1) ? *a & FILLMASK : 1;
			}
			incr_a = false;
		}
		if (incr_b) {
			while(b_pos <= res_pos) {
				++b;
				b_pos += (*b & BIT1) ? *b & FILLMASK : 1;
			}
			incr_b = false;
		}
		if (a_pos == b_pos) {
			if ((*a & ONEFILL) == ONEFILL && (*b & ONEFILL) == ONEFILL)
				next_word = ONEFILL | (a_pos - res_pos);
			else if ((*a & BIT1) || (*b & BIT1))
				next_word = BIT1 | (a_pos - res_pos);
			else {
				uint64_t u = *a & *b;
				next_word = (u == 0) ? BIT1 | 1 : u;
			}
			incr_a = true;
			incr_b = true;
			res_pos = a_pos;
		}
		else if (a_pos < b_pos) {
			if ((*b & ONEFILL) == ONEFILL) {
				if((*a & ONEFILL) == ONEFILL)
					next_word = ONEFILL | (a_pos - res_pos);
				else if (*a & BIT1)
					next_word = BIT1 | (a_pos - res_pos);
				else
					next_word = *a;
				res_pos = a_pos;
				incr_a = true;
			}
			else { // b is a 0-fill because it can't be a literal word and have b_pos > a_pos
				next_word = BIT1 | (b_pos - res_pos);
				res_pos = b_pos;
				incr_a = true;
				incr_b = true;
			}
		}
		else { // a_pos > b_pos
			if ((*a & ONEFILL) == ONEFILL) {
				if((*b & ONEFILL) == ONEFILL)
					next_word = ONEFILL | (b_pos - res_pos);
				else if (*b & BIT1)
					next_word = BIT1 | (b_pos - res_pos);
				else
					next_word = *b;
				res_pos = b_pos;
				incr_b = true;
			}
			else { // a is a 0-fill because it can't be a literal word and have a_pos > b_pos
				next_word = BIT1 | (a_pos - res_pos);
				res_pos = a_pos;
				incr_a = true;
				incr_b = true;
			}
		}
		if ((next_word & BIT1)
			&& res.size() > 0
			&& (res.back() & BIT1)
			&& (res.back() & BIT2) == (next_word & BIT2)
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

