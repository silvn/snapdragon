#include "bvec.h"

// constructor - given a sorted vector of distinct 32bit integers
BitVector::BitVector(vector<word_t>& vals) {
    count = vals.size();
    // if the density is too low, run length encoding will take MORE space
    if (lowDensity(vals)) {
        if (DEBUG) printf("constructed non rle\n");
        words = vals;
        rle = false;
        size = words.back();
    }
    else {
        if (DEBUG) printf("constructRLE\n");
        constructRLE(vals);
    }
}

vector<word_t>& BitVector::getWords() {
    return words;
}

word_t BitVector::getSize() { return size; }
word_t BitVector::bytes() { return 4 * words.size(); }

void BitVector::compress() {
    if (rle) { /* Throw exception? */ return; }
    vector<word_t> tmp;
    tmp.swap(words);
    constructRLE(tmp);
}

// in place version of the bitwise OR operator.
void BitVector::operator|=(BitVector& bv) {
    // decide which version we'll be using
    if (rle)
        if (bv.rle)
            rleORrle(bv);
        else
            rleORnon(bv);
    else
        if (bv.rle)
            nonORrle(bv);
        else
            nonORnon(bv);
}
// in place version of the bitwise AND operator.
void BitVector::operator&=(BitVector& bv) {
    // decide which version we'll be using
    if (rle)
        if (bv.rle)
            rleANDrle(bv);
        else
            rleANDnon(bv);
    else
        if (bv.rle)
            nonANDrle(bv);
        else
            nonANDnon(bv);
}
BitVector* BitVector::operator|(BitVector& rhs) {
    BitVector *res = new BitVector();
    if (words.size() > rhs.words.size()) {
        res->copy(rhs);
        *res |= *this;
        return res;
    }
    res->copy(*this);
    *res |= rhs;
    return res;
}
BitVector* BitVector::operator&(BitVector& rhs) {
    BitVector *res = new BitVector();
    if (words.size() > rhs.words.size()) {
        res->copy(rhs);
        *res &= *this;
        return res;
    }
    res->copy(*this);
    *res &= rhs;
    return res;
}

bool BitVector::operator==(BitVector& other) const {
    return (words == other.words) &&
           (count == other.count) &&
           (size  == other.size)  &&
           (rle   == other.rle);
}

bool BitVector::equals(const BitVector& other) const {
    return (words == other.words) &&
           (count == other.count) &&
           (size  == other.size)  &&
           (rle   == other.rle);
}

BitVector* BitVector::copyflip() {
    BitVector *res = new BitVector();
    res->copy(*this);
    res->flip();
    return res;
}

void BitVector::nonORnon(BitVector& bv) {
    
    vector<word_t> res;
    vector<word_t>::iterator a = words.begin();
    vector<word_t>::iterator b = bv.words.begin();
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

void BitVector::nonANDnon(BitVector& bv) {
    
    vector<word_t> res;
    vector<word_t>::iterator a = words.begin();
    vector<word_t>::iterator b = bv.words.begin();
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

void BitVector::nonANDrle(BitVector& bv) {
    // decompress
    // run nonANDnon
    BitVector *tmp = new BitVector();
    tmp->copy(bv);
    tmp->decompress();
    nonANDnon(*tmp);
    delete tmp;
}

void BitVector::rleANDnon(BitVector& bv) {
    // decompress
    // run nonANDnon
    decompress();
    nonANDnon(bv);
}

void BitVector::nonORrle(BitVector& bv) {
    // compress
    // run rleORrle
    compress();
    rleORrle(bv);
}

void BitVector::rleORnon(BitVector& bv) {
    // compress
    // run rleORrle
    BitVector *tmp = new BitVector();
    tmp->copy(bv);
    tmp->compress();
    rleORrle(*tmp);
    delete tmp;
}

void BitVector::setBit(word_t x) {
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
            vector<word_t>::iterator lb = lower_bound(words.begin(),words.end(),x);
            if (*lb == x) return;
            words.insert(lb,x);
        }
        size++;
        count++;
    }
}

BitVector&
BitVector::copy(const BitVector& bv) {
    words = bv.words;
    count = bv.count;
    size = bv.size;
    rle = bv.rle;
    return *this;
}

ostream& operator<<(ostream &os, const vector<word_t> vec) {
    return os << vec;
}

ostream& operator<<(ostream &os, const BitVector &vec) {
    return os << vec;
}

void
BitVector::save(const BitVector &bv, const char *filename) {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << bv;
}

void
BitVector::restore(BitVector &bv, const char *filename) {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> bv;
}
