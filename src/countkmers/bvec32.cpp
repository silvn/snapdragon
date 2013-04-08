#include "bvec32.h"

// constructor - given a sorted vector of distinct 32bit integers
bvec32::bvec32(vector<uint32_t>& vals) {
	count = vals.size();
	// if the density is too low, run length encoding will take MORE space
	if (low_density(vals)) {
		if (DEBUG) printf("constructed non rle\n");
		words = vals;
		rle = false;
		size = words.back();
	}
	else {
		if (DEBUG) printf("construct_rle\n");
		construct_rle(vals);
	}
}

// constructor - given a previously dumped bvec
bvec32::bvec32(uint32_t *buf) {
	uint32_t nwords = buf[0];
	size = buf[1];
	count = buf[2];
	rle = false;
	if (nwords & BIT1) {
		nwords -= BIT1;
		rle=true;
	}
	words.resize(nwords);
	memcpy(words.data(),buf+3,nwords*4);
}

vector<uint32_t>& bvec32::get_words() {
	return words;
}
// DIY serialization
size_t bvec32::dump(uint32_t **buf) {
	// allocate space in buf
	size_t dbytes = sizeof(uint32_t)*(3 + words.size());
	*buf = (uint32_t*)malloc(dbytes);
	if (*buf == NULL) {
		fprintf(stderr,"failed to allocate %zi bytes\n",dbytes);
		return 0;
	}
	(*buf)[0] = words.size();
	if (rle) (*buf)[0] |= BIT1;
	(*buf)[1] = size;
	(*buf)[2] = count;
	memcpy(*buf + 3, words.data(), bytes());
	return dbytes;
}

bool bvec32::low_density(vector<uint32_t>& vals) {
	if (DEBUG) printf("low_density() %u/%u %c %f\n", (uint32_t)vals.size(),
		vals.back() - vals.front() + 1,
			((double)vals.size()/(double)(vals.back() - vals.front() + 1) < 1.0/(double)LITERAL_SIZE) ? '<' : '>',
				1.0/(double)LITERAL_SIZE);
	return (double)vals.size()/(double)(vals.back() - vals.front() + 1) < 1.0/(double)LITERAL_SIZE;
	
}

void bvec32::print() {
	printf("rle: %c\n",rle ? 'T' : 'F');
	printf("words:\n");
	for(int i=0;i<words.size();i++) {
		printf(" %i ",i);
		if (rle)
			if ((words[i] & ONEFILL) == ONEFILL) 
				printf("1-fill %zi %zi\n",words[i] & FILLMASK, LITERAL_SIZE*(words[i] & FILLMASK));
			else if (words[i] & BIT1)
				printf("0-fill %zi %zi\n",words[i] & FILLMASK, LITERAL_SIZE*(words[i] & FILLMASK));
			else
				printf("literal %u\n",words[i]);
		else
			printf("direct %u\n",words[i]);
	}
}

void bvec32::construct_rle(vector<uint32_t>& vals) {
	rle = true;
	if (vals.size() == 0) {
		size=0;
		count=0;
		words.clear();
		words.push_back(0); // current word is an empty literal
		return;
	}
	uint32_t word_end = LITERAL_SIZE - 1;
	uint32_t word=0;
	uint32_t gap_words = vals.front()/LITERAL_SIZE;
	if (gap_words > 0) {
		word_end += LITERAL_SIZE*gap_words;
		while (gap_words > FILLMASK) {
			words.push_back(ZEROFULL);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			words.push_back(gap_words | BIT1);
	}
	for(vector<uint32_t>::iterator ii = vals.begin(); ii != vals.end(); ++ii) {
		if (*ii == word_end)
			word |= 1;
		else if (*ii < word_end)
			word |= ((uint32_t)1 << (word_end - *ii));
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
				? 1 : (uint32_t)1 << (word_end - *ii);
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

void bvec32::compress() {
	if (rle) { /* Throw exception? */ return; }
	vector<uint32_t> tmp;
	tmp.swap(words);
	construct_rle(tmp);
	rle = true;
}

void bvec32::decompress() {
	if (!rle) { /* Throw exception? */ return; }
	// retrieve the set bits from the compressed vector
	vector<uint32_t> res;
	res.reserve(cnt());
	uint32_t pos=0;
	for(vector<uint32_t>::iterator ii = words.begin(); ii != words.end(); ++ii) {
		if ((*ii & ONEFILL) == ONEFILL) {
			uint32_t n_ones = LITERAL_SIZE*(*ii & FILLMASK);
			for(uint32_t i=0;i < n_ones; i++) {
				res.push_back(pos);
				pos++;
			}
		}
		else if (*ii & BIT1)
			pos += LITERAL_SIZE*(*ii & FILLMASK);
		else { // desconstruct the literal word
			for(uint32_t bit=1; bit<=LITERAL_SIZE; bit++)
				if (*ii & ((uint32_t)1 << (LITERAL_SIZE-bit)))
					res.push_back(pos+bit-1);
			pos += LITERAL_SIZE;
		}
	}
	words.swap(res);
	count = words.size();
	rle = false;
}

// in place version of the bitwise OR operator.
void bvec32::operator|=(bvec32& bv) {
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
void bvec32::operator&=(bvec32& bv) {
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
bvec32* bvec32::operator|(bvec32& rhs) {
	bvec32 *res = new bvec32();
	if (words.size() > rhs.words.size()) {
		res->copy(rhs);
		*res |= *this;
		return res;
	}
	res->copy(*this);
	*res |= rhs;
	return res;
}
bvec32* bvec32::operator&(bvec32& rhs) {
	bvec32 *res = new bvec32();
	if (words.size() > rhs.words.size()) {
		res->copy(rhs);
		*res &= *this;
		return res;
	}
	res->copy(*this);
	*res &= rhs;
	return res;
}

bool bvec32::operator==(bvec32& other) const {
	return (words == other.words) &&
		   (count == other.count) &&
		   (size  == other.size)  &&
		   (rle	  == other.rle);
}

bool bvec32::equals(const bvec32& other) const {
	return (words == other.words) &&
		   (count == other.count) &&
		   (size  == other.size)  &&
		   (rle	  == other.rle);
}

void bvec32::non_OR_non(bvec32& bv) {
	
	vector<uint32_t> res;
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();
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

void bvec32::non_AND_non(bvec32& bv) {
	
	vector<uint32_t> res;
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();
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

void bvec32::non_AND_rle(bvec32& bv) {
	// decompress
	// run non_AND_non
	bvec32 *tmp = new bvec32();
	tmp->copy(bv);
	tmp->decompress();
	non_AND_non(*tmp);
	delete tmp;
}
void bvec32::rle_AND_non(bvec32& bv) {
	// decompress
	// run non_AND_non
	decompress();
	non_AND_non(bv);
}
void bvec32::non_OR_rle(bvec32& bv) {
	// compress
	// run rle_OR_rle
	compress();
	rle_OR_rle(bv);
}
void bvec32::rle_OR_non(bvec32& bv) {
	// compress
	// run rle_OR_rle
	bvec32 *tmp = new bvec32();
	tmp->copy(bv);
	tmp->compress();
	rle_OR_rle(*tmp);
	delete tmp;
}

void bvec32::matchSize(bvec32 &bv) {
	if (size < bv.size) {
		uint32_t gap_words = (bv.size - size)/LITERAL_SIZE;
		while (gap_words > FILLMASK) {
			words.push_back(ZEROFULL);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			words.push_back(BIT1 | gap_words);
		size = bv.size;
	}
	else if (size > bv.size) {
		uint32_t gap_words = (size - bv.size)/LITERAL_SIZE;
		while (gap_words > FILLMASK) {
			bv.words.push_back(ZEROFULL);
			gap_words -= FILLMASK;
		}
		if (gap_words > 0)
			bv.words.push_back(BIT1 | gap_words);
		bv.size = size;
	}
}

void bvec32::rle_OR_rle(bvec32& bv) {
	// ensure that both bvec32s are the same size
	this->matchSize(bv);
	if (size == 0)
		return;
	
	vector<uint32_t> res; // fill this then swap with this.words
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();

	// maintain the end position of the current word
	uint32_t a_pos = (*a & BIT1) ? (*a & FILLMASK) : 1;
	uint32_t b_pos = (*b & BIT1) ? (*b & FILLMASK) : 1;
	uint32_t res_pos=0;
	uint32_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	uint32_t last_pos = size/LITERAL_SIZE;
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
			uint32_t n_words = (res.back() & FILLMASK) + (next_word & FILLMASK);
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
void bvec32::rle_AND_rle(bvec32& bv) {
	// ensure that both bvec32s are the same size
	this->matchSize(bv);
	if (size == 0)
		return;
	
	vector<uint32_t> res; // fill this then swap with this.words
	vector<uint32_t>::iterator a = words.begin();
	vector<uint32_t>::iterator b = bv.words.begin();

	// maintain the end position of the current word
	uint32_t a_pos = (*a & BIT1) ? *a & FILLMASK : 1;
	uint32_t b_pos = (*b & BIT1) ? *b & FILLMASK : 1;
	uint32_t res_pos=0;
	uint32_t next_word;
	bool incr_a = false;
	bool incr_b = false;
	uint32_t last_pos = size/LITERAL_SIZE;
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
				uint32_t u = *a & *b;
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
			uint32_t n_words = (res.back() & FILLMASK) + (next_word & FILLMASK);
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

bool bvec32::find(uint32_t x) {
	if (rle)
		return rle_find(x);
	else
		return binary_search(words.begin(), words.end(), x);
}

bool bvec32::rle_find(uint32_t x) {
	uint32_t pos=0;
	for(vector<uint32_t>::iterator it=words.begin();it!=words.end();++it) {
		// what type of word is it?
		if (*it & BIT1) { // fill word
			pos += (*it & FILLMASK) * LITERAL_SIZE;
			if (pos >= x) {
				if (*it & BIT2) // 1-fill
					return true;
				return false;
			}
		}
		else { // literal word
			pos += LITERAL_SIZE;
			if (pos >= x) {
				if ((1UL << (pos-x)) & *it)
					return true;
				return false;
			}
		}
	}
	return false;
}

void bvec32::setBit(uint32_t x) {
	if (rle) { // this should really be done in an rle specific way.
		decompress();
		setBit(x);
		compress();
	}
	else {
		if (words.size() == 0 || x > words.back()) {
			words.push_back(x);
		}
		else {
			vector<uint32_t>::iterator lb = lower_bound(words.begin(),words.end(),x);
			if (*lb == x) return;
			words.insert(lb,x);
		}
		size++;
		count++;
	}
}

// size holds the last meaningful bit although the last word may contain space for other bits
void bvec32::appendFill(bool bit, uint32_t n) {
	if (!rle) compress();
//	fprintf(stderr,"appendFill(%i,%u) words.size(): %zi\n",bit ? 1 : 0 , n, words.size());
	if (n==0) return;
//	print();
	if ((words.back() & BIT1) == 0) { // current word is a literal word
		// size % LITERAL_SIZE is the number of bits already used
		uint32_t bits_available = LITERAL_SIZE - size % LITERAL_SIZE;
//		fprintf(stderr,"bits_available: %u\n",bits_available);
		if(bit) {
			// append some ones
			uint32_t append = (1 << bits_available) - 1;
			if (n < bits_available) {
				append = ((1 << n) - 1) << (bits_available - n);
				size += n;
				n=0;
			}
			else {
				size += bits_available;
				n -= bits_available;
				// the word is full, check if we should convert to a 1-fill
				if (words.back() == ALL1S) {
					words.back() = ONEFILL1;
				}
			}
			words.back() |= append;
		}
		else {
			if (n < bits_available) {
				size += n;
				n=0;
			}
			else {
				size += bits_available;
				n -= bits_available;
				// the word is full, check if we should convert to a 0-fill
				if (words.back() == 0) {
					words.back() = ZEROFILL1;
				}
			}
		}
	}
//	print();
	if (n==0) return;
	// append/update fill words
	uint32_t n_fills = n/LITERAL_SIZE;
	if (n_fills > 0) {
		if (bit) {
			if (words.back() & ONEFILL) { // extend previous 1-fill
				words.back() += n_fills;
			}
			else { // append a 1-fill
				words.push_back(ONEFILL | n_fills);
			}
		}
		else {
			if (words.back() & BIT1 && !(words.back() & BIT2)) { // extend previous 0-fill
				words.back() += n_fills;
			}
			else { // append a 1-fill
				words.push_back(BIT1 | n_fills);
			}
		}
		n -= n_fills*LITERAL_SIZE;
		size += n_fills*LITERAL_SIZE;
	}
//	print();
	if (n==0) return;
	// add the remaining bits to a literal word
	if (bit)
		words.push_back(((1<<n)-1) << (LITERAL_SIZE - n));
	else
		words.push_back(0);
	size += n;
//	print();
}
