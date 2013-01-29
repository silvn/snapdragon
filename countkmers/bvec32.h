#ifndef BVEC32_H
#define BVEC32_H

#define WORD_SIZE 32
#define LITERAL_SIZE 31
#define BIT1 0x80000000
#define BIT2 0x40000000
#define FILLMASK 0x3FFFFFFF
#define ALL1S 0x7FFFFFFF
#define ONEFILL 0xC0000000
#define ONEFILL1 0xC0000001
#define ONEFULL 0xFFFFFFFF
#define ZEROFULL 0xAFFFFFFF

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
	bvec(vector<uint32_t>& vals);
	void print();
	void compress();
    void decompress();

	// logical set operations
	void operator|=(bvec& rhs);
	bvec* operator|(bvec&);
	void operator&=(bvec& rhs);
	bvec* operator&(bvec&);

	// is x in the set?
	bool find(uint32_t x);

	// basic metrics
	inline uint32_t cnt();
	inline uint32_t bytes() { return 4*words.size(); }

private:
	vector<uint32_t> words;
	bool rle;
	uint32_t count; // cache the number of set bits
	uint32_t size; // bits in the uncompressed bitvector

	bool low_density(vector<uint32_t>& vals);
	void construct_rle(vector<uint32_t>& vals);
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
inline uint32_t bvec::cnt() {
	count=0;
	if (count == 0) 
		if (rle)
			for(vector<uint32_t>::iterator it = words.begin(); it != words.end(); ++it)
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