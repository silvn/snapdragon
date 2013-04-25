#ifndef SNAPDRAGON_BVEC_H
#define SNAPDRAGON_BVEC_H

#define DEBUG false

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <boost/serialization/vector.hpp>
#include <boost/cstdint.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace std;

typedef uint32_t word_t;

class BitVector {
    vector<word_t> words;
    bool rle;
    word_t count; // cache the number of set bits
    word_t size; // bits in the uncompressed bitvector

    friend class boost::serialization::access;
    friend ostream & operator <<(ostream &, const BitVector &);

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & words & count & size & rle;
    }

    // checkpoint: active word and associated info
    struct checkpoint {
        vector<word_t>::iterator active_word;
        word_t bit_pos;
    } frontier;

public:
    // Destructor
    ~BitVector() {};

    // Constructors
    BitVector()                     : rle(false), count(0), size(0) {};
    BitVector(bool wah)             : rle(false), count(0), size(0) {
        compress();
    };
    BitVector(vector<word_t>& vals);
    BitVector(word_t* buf); // DIY deserialization
    
    BitVector& copy(const BitVector& bv);
    
    void   print();
    void   compress();
    void   decompress();
    size_t dump(word_t **buf); // DIY serialization

    // logical set operations
    void                 flip();
    BitVector *          copyflip();
    void        operator |= (BitVector& rhs);
    BitVector * operator |  (BitVector&);
    void        operator &= (BitVector& rhs);
    BitVector * operator &  (BitVector&);
    bool        operator ==(BitVector&) const;
    
    bool                 equals(const BitVector&) const;
    vector<word_t>&      getWords();

    // is x in the set?
    bool find(word_t x);

    // find the position of the next set bit after x
    word_t nextOne(word_t x);

    // insert x into an existing BitVector (at the end is faster)
    void setBit(word_t x);

    // for constructing a rle BitVector one bit at a time
    void appendFill(bool bit, word_t count);

    // basic metrics
    word_t cnt();
    word_t getSize();
    word_t bytes();
    
    static void save(const BitVector &bv, const char *filename);
    static void restore(BitVector &bv, const char *filename);

private:

    bool   lowDensity(vector<word_t>& vals);
    void   constructRLE(vector<word_t>& vals);
    void   matchSize(BitVector& bv);
    void   rleORrle(BitVector& rhs);
    void   rleORnon(BitVector& rhs);
    void   nonORrle(BitVector& rhs);
    void   nonORnon(BitVector& rhs);
    void   rleANDrle(BitVector& rhs);
    void   rleANDnon(BitVector& rhs);
    void   nonANDrle(BitVector& rhs);
    void   nonANDnon(BitVector& rhs);
    word_t popCount(word_t val) const;

};

#endif
