#ifndef BVEC_H
#define BVEC_H
#include <algorithm>
#include <vector>
using namespace std;

#define WORD_SIZE 32
#define LITERAL_SIZE 31

class bvec
{
public:
	// Desctructor
	~bvec() {};

	// Constructors
	bvec();
	bvec(vector<uint32_t>& vals);
	bvec(vector<uint64_t>& vals);

	// logical set operations
	bvec* operator&(const bvec& rhs); // AND
	bvec* operator|(const bvec& rhs); // OR
	bvec* operator^(const bvec& rhs); // XOR
	bvec* operator-(bvec& rhs); // subtract
	void operator&=(const bvec& rhs);
	void operator|=(const bvec& rhs);
	void operator^=(const bvec& rhs);
	void operator-=(bvec& rhs);
	void flip(); // NOT

	// count the number of values stored
	inline uint64_t cnt();
	// is x in the set?
	bool find(uint64_t x);

private:
	vector<uint32_t> words;
	uint64_t count; // cache the number of set bits

	uint32_t litbit  = 1 << 31; // literal bit: 0 is a fill word, 1 is a literal word
	uint32_t fillbit = 1 << 30; // fill bit is 1 if the fill word is all ones
	uint32_t fillmask = fillbit - 1; // these bits say how long the fill is in 31mers

	inline uint64_t fill_length(uint32_t word);
};


inline uint32_t bvec::cnt() {
	if (count == 0) 
		for(vector<uint32_t>::iterator it = words.begin(); it != words.end(); ++it)
			count += (*it & litbit) ? __builtin_popcount(*it) - 1 : (*it & fillbit) ? (*it & fillmask) * LITERAL_SIZE  : 0;
	return count;
}


#endif