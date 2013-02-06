#ifndef BVEC32_H
#define BVEC32_H

#define DEBUG false

#define WORD_SIZE 32
#define LITERAL_SIZE 31
#define BIT1     0x80000000UL
#define BIT2     0x40000000UL
#define FILLMASK 0x3FFFFFFFUL
#define ALL1S    0x7FFFFFFFUL
#define ONEFILL  0xC0000000UL
#define ONEFILL1 0xC0000001UL
#define ONEFULL  0xFFFFFFFFUL
#define ZEROFULL 0xAFFFFFFFUL

#include <algorithm>
#include <vector>
#include <boost/serialization/vector.hpp>
using namespace std;

class bvec32 {
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & rle & count & size & words;
	}
	
public:
	// Destructor
	~bvec32() {};

	// Constructors
	bvec32() {};
	bvec32(vector<uint32_t>& vals);
	void print();
	void compress();
    void decompress();
    vector<uint32_t>& get_words();

	// logical set operations
	void operator|=(bvec32& rhs);
	bvec32* operator|(bvec32&);
	void operator&=(bvec32& rhs);
	bvec32* operator&(bvec32&);

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
	inline bvec32& copy(const bvec32& bv);
	void matchSize(bvec32& bv);
	void rle_OR_rle(bvec32& rhs);
	void rle_OR_non(bvec32& rhs);
	void non_OR_rle(bvec32& rhs);
	void non_OR_non(bvec32& rhs);
	void rle_AND_rle(bvec32& rhs);
	void rle_AND_non(bvec32& rhs);
	void non_AND_rle(bvec32& rhs);
	void non_AND_non(bvec32& rhs);
	inline uint32_t popcount(uint32_t val) const;
};

inline uint32_t bvec32::popcount(uint32_t val) const {
	// number of 1 bits in a value between 0 and 255
	static const uint32_t table[256] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
    return table[val&0xFFUL] + table[(val>>8)&0xFFUL] +
	table[(val>>16)&0xFFUL] + table[(val>>24)&0xFFUL];
}


// count the number of set bits
inline uint32_t bvec32::cnt() {
	if (count == 0)
		if (rle)
			for(vector<uint32_t>::iterator it = words.begin(); it != words.end(); ++it)
				count += (*it & BIT1) ? (*it & BIT2) ? (*it & FILLMASK) * LITERAL_SIZE : 0 : popcount(*it & ALL1S);
		else
			count = words.size();
	return count;
}

// count the number of set bits
// inline uint32_t bvec32::bvcnt() {
// 	count=0;
// 	if (count == 0) 
// 		if (rle)
// 			for(vector<uint32_t>::iterator it = words.begin(); it != words.end(); ++it)
// 				count += (*it & BIT1) ? (*it & BIT2) ? (*it & FILLMASK) * LITERAL_SIZE : 0 : __builtin_popcount(*it);
// 		else
// 			count = words.size();
// 	return count;
// }

// make a copy
inline bvec32& bvec32::copy(const bvec32& bv) {
	words = bv.words;
	count = bv.count;
	size = bv.size;
	rle = bv.rle;
	return *this;
}

#endif