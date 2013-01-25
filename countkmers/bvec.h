#ifndef BVEC_H
#define BVEC_H
#include <algorithm>
#include <vector>
//#include <boost/serialization/vector.hpp>
using namespace std;

#define WORD_SIZE 32
#define LITERAL_SIZE 31
#define BIT32 2147483648
#define BIT31 1073741824
#define FILLMASK 1073741823
#define ALL1S 2147483647
#define ONEFILL 3221225472
#define ONEFILL1 3221225473
#define ONEFILL2 4294967295

class bvec
{
	// friend class boost::serialization::access;

	// template<class Archive>
	// void serialize(Archive & ar, const unsigned int version) {
	// 	ar & words;
	// }
	
public:
	// Destructor
	~bvec() {};

	// Constructors
	bvec();
	bvec(vector<uint32_t>& vals);
	bvec(vector<uint64_t>& vals);

	// logical set operations
	bvec& union(bvec& bv);
	bvec& intersect(bvec& bv);

	// is x in the set?
	bool find(uint64_t x);

	// basic metrics
	inline uint64_t cnt();
	inline uint32_t bytes() { return 4*words.size(); }

private:
	vector<uint32_t> words;
	uint64_t count; // cache the number of set bits
	
};

// count the number of set bits
inline uint64_t bvec::cnt() {
	if (count == 0) 
		for(vector<uint32_t>::iterator it = words.begin(); it != words.end(); ++it)
			count += (*it & BIT32) ? (*it & BIT31) ? (*it & FILLMASK) * LITERAL_SIZE : 0 : __builtin_popcount(*it) - 1;
	return count;
}


#endif