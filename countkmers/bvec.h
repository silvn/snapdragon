#ifndef BVEC_H
#define BVEC_H
#include <algorithm>
#include <vector>
#include <boost/serialization/vector.hpp>
using namespace std;

#define WORD_SIZE 32
#define LITERAL_SIZE 31
#define BIT32 2147483648
#define BIT31 1073741824
#define FILLMASK 1073741823
#define ALL1S 2147483647
#define ONEFILL 3221225472
#define ONEFILL1 3221225473

class bvec
{
	friend class boost::serialization::access;
	vector<uint32_t> words;
	uint64_t count; // cache the number of set bits

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & words;
	}
	
public:
	// Desctructor
	~bvec() {};

	// Constructors
	bvec();
	bvec(vector<uint32_t>& vals);
	bvec(vector<uint64_t>& vals);

	// basic metrics
	uint64_t cnt();
	uint32_t bytes() { return 4*words.size(); }
	uint32_t size() { } // uncompressed length

	// logical set operations

	// is x in the set?
	bool find(uint64_t x);
};


#endif