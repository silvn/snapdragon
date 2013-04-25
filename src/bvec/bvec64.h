#ifndef BVEC64_H
#define BVEC64_H

#define DEBUG false

#define WORD_SIZE 64
#define LITERAL_SIZE 63
#define BIT1     0x8000000000000000
#define BIT2     0x4000000000000000
#define FILLMASK 0x3FFFFFFFFFFFFFFF
#define ALL1S    0x7FFFFFFFFFFFFFFF
#define ONEFILL  0xC000000000000000
#define ONEFILL1 0xC000000000000001
#define ONEFULL  0xFFFFFFFFFFFFFFFF
#define ZEROFULL 0xAFFFFFFFFFFFFFFF

#include <algorithm>
#include <vector>
#include <boost/serialization/vector.hpp>
using namespace std;

class bvec64 {
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & rle & count & size & words;
	}
	
public:
	// Destructor
	~bvec64() {};

	// Constructors
	bvec64() {};
	bvec64(vector<uint64_t>& vals);
	void print();
	void compress();
    void decompress();
    vector<uint64_t>& getWords();

	// logical set operations
	void operator|=(bvec64& rhs);
	bvec64* operator|(bvec64&);
	void operator&=(bvec64& rhs);
	bvec64* operator&(bvec64&);

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

	bool lowDensity(vector<uint64_t>& vals);
	void constructRLE(vector<uint64_t>& vals);
	inline bvec64& copy(const bvec64& bv);
	void matchSize(bvec64& bv);
	void rleORrle(bvec64& rhs);
	void rleORnon(bvec64& rhs);
	void nonORrle(bvec64& rhs);
	void nonORnon(bvec64& rhs);
	void rleANDrle(bvec64& rhs);
	void rleANDnon(bvec64& rhs);
	void nonANDrle(bvec64& rhs);
	void nonANDnon(bvec64& rhs);

	inline uint64_t popCount(uint64_t val) const;
};

inline uint64_t bvec64::popCount(uint64_t val) const {
	// number of 1 bits in a value between 0 and 255
	static const uint64_t table[256] = {
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
    return table[val&0xFFULL] + table[(val>>8)&0xFFULL] +
	table[(val>>16)&0xFFULL] + table[(val>>24)&0xFFULL] +
	table[(val>>32)&0xFFULL] + table[(val>>40)&0xFFULL] +
	table[(val>>48)&0xFFULL] + table[(val>>56)&0xFFULL];
}

// count the number of set bits
inline uint64_t bvec64::cnt() {
	if (count == 0)
		if (rle)
			for(vector<uint64_t>::iterator it = words.begin(); it != words.end(); ++it)
				count += (*it & BIT1) ? (*it & BIT2) ? (*it & FILLMASK) * LITERAL_SIZE : 0 : popCount(*it & ALL1S);
		else
			count = words.size();
	return count;
}

// make a copy
inline bvec64& bvec64::copy(const bvec64& bv) {
	words = bv.words;
	count = bv.count;
	size = bv.size;
	rle = bv.rle;
	return *this;
}

#endif
