#ifndef BVEC64_H
#define BVEC64_H

#define DEBUG false

#define WORD_SIZE 64
#define LITERAL_SIZE 63
#define BIT1 0x8000000000000000
#define BIT2 0x4000000000000000
#define FILLMASK 0x3FFFFFFFFFFFFFFF
#define ALL1S 0x7FFFFFFFFFFFFFFF
#define ONEFILL 0xC000000000000000
#define ONEFILL1 0xC000000000000001
#define ONEFULL 0xFFFFFFFFFFFFFFFF
#define ZEROFULL 0xAFFFFFFFFFFFFFFF

#include <algorithm>
#include <vector>
#include <boost/serialization/vector.hpp>
using namespace std;

class bvec {
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & rle & count & size & words;
	}
	
public:
	// Destructor
	~bvec() {};

	// Constructors
	bvec() {};
	bvec(vector<uint64_t>& vals);
	void print();
	void compress();
    void decompress();
    vector<uint64_t>& get_words();

	// logical set operations
	void operator|=(bvec& rhs);
	bvec* operator|(bvec&);
	void operator&=(bvec& rhs);
	bvec* operator&(bvec&);

	// is x in the set?
	bool find(uint64_t x);

	// basic metrics
	inline uint64_t cnt();
	inline uint64_t bytes() { return 8*words.size(); }

private:
	vector<uint64_t> words;
	bool rle;
	uint64_t count; // cache the number of set bits
	uint64_t size; // bits in the uncompressed bitvector

	bool low_density(vector<uint64_t>& vals);
	void construct_rle(vector<uint64_t>& vals);
	inline bvec& copy(const bvec& bv);
	void matchSize(bvec& bv);
	void rle_OR_rle(bvec& rhs);
	void rle_OR_non(bvec& rhs);
	void non_OR_rle(bvec& rhs);
	void non_OR_non(bvec& rhs);
	void rle_AND_rle(bvec& rhs);
	void rle_AND_non(bvec& rhs);
	void non_AND_rle(bvec& rhs);
	void non_AND_non(bvec& rhs);
};

// count the number of set bits
inline uint64_t bvec::cnt() {
	count=0;
	if (count == 0) 
		if (rle)
			for(vector<uint64_t>::iterator it = words.begin(); it != words.end(); ++it)
				count += (*it & BIT1) ? (*it & BIT2) ? (*it & FILLMASK) * LITERAL_SIZE : 0 : __builtin_popcount(*it);
		else
			count = words.size();
	return count;
}

// make a copy
inline bvec& bvec::copy(const bvec& bv) {
	words = bv.words;
	count = bv.count;
	size = bv.size;
	rle = bv.rle;
	return *this;
}

#endif
