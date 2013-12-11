#ifndef SNAPDRAGON_BITVECTOR_H
#define SNAPDRAGON_BITVECTOR_H
#include <cstdio>
#include <vector>
using namespace std;
#define ZEROFILL 0
#define ONEFILL 1
#define LITERAL 2
#define NONEOFTHEABOVE 3

template <class T>
class BitVector {
public:
    BitVector(); // constructor
    BitVector(T *buf); // reconstructor
    size_t dump(T **buf); // serialize

    void appendWord(T word); // append uncompressed bits in word 
    void appendFill0(size_t length); // extend the bitvector by length 0 bits
    void appendFill1(size_t length); // extend the bitvector by length 1 bits
    void inflateWord(T *word, size_t wordStart); // fills a word with uncompressed bits starting at wordStart
    
    uint64_t getSize() { return size; }
    uint64_t getCount() { return count; }

    BitVector<T>* operator&(BitVector<T>* rhs);
    BitVector<T>* operator|(BitVector<T>* rhs);
    void operator&=(BitVector<T>* rhs);
    void operator|=(BitVector<T>* rhs);
    void flip();

private:
    vector<T> words;  // a mix of literal and fill words. MSB of fill words indicate the type of fill
    vector<T> isFill; // uncompressed bitvector indicating which words are fills

    uint64_t count;   // cache the number of set bits
    uint64_t size;    // bits in the uncompressed bitvector
    int nbits;        // number of bits per LITERAL word sizeof(T)*8
    int shiftby;      // log base 2 of nbits
    T onefill1;       // e.g., 1000000000000000001

    // active word used for iterating
    
    size_t   activeWordIdx;   // offset in to words vector
    char     activeWordType;  // {ONEFILL, ZEROFILL, LITERAL, NONEOFTHEABOVE}
    uint64_t activeWordStart; // uncompressed bit position at start of word
    uint64_t activeWordEnd;   // uncompressed bit position after last bit in a word

    void seek(size_t wordStart); // locate the activeWord that contains wordStart
    void firstActiveWord(); // jump directly to first word
    void nextActiveWord(); // advance to next word
};

// template <class T>
// void LCBitSlicedIndex<T>::transpose(uint64_t A[64]) {
void transpose(uint64_t A[64]) {
    int j, k;
    uint64_t m, t;
    m = 0x00000000FFFFFFFFULL;
    for (j = 32; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 64; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}
// template <class T>
// void LCBitSlicedIndex<T>::transpose(uint32_t A[32]) {
void transpose(uint32_t A[32]) {
    int j, k;
    uint32_t m, t;
    
    m = 0x0000FFFF;
    for (j = 16; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 32; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}
// template <class T>
// void LCBitSlicedIndex<T>::transpose(uint16_t A[16]) {
void transpose(uint16_t A[16]) {
    int j, k;
    uint16_t m, t;
    
    m = 0x00FF;
    for (j = 8; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 16; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}

// for debugging
void printBits(size_t const size, void const * const ptr) {
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    for (i=size-1;i>=0;i--)
        for (j=7;j>=0;j--){
            byte = b[i] & (1<<j);
            byte >>= j;
           fprintf(stderr,"%u", byte);
        }
}

template <class T>
BitVector<T>::BitVector() {
    nbits = sizeof(T)*8;
    onefill1 = ((T)1 << (nbits-1)) | (T)1;
    if      (nbits == 64) shiftby = 6;
    else if (nbits == 32) shiftby = 5;
    else if (nbits == 16) shiftby = 4;
    else shiftby = 3;
    count = 0;
    size  = 0;
    activeWordIdx   = 0;
    activeWordStart = 0;
    activeWordEnd   = 0;
    activeWordType  = NONEOFTHEABOVE;
}

// constructor - given a previously dumped BitVector
template <class T>
BitVector<T>::BitVector(T *buf) {
    nbits = sizeof(T)*8;
    onefill1 = ((T)1 << (nbits-1)) | (T)1;
    if      (nbits == 64) shiftby = 6;
    else if (nbits == 32) shiftby = 5;
    else if (nbits == 16) shiftby = 4;
    else shiftby = 3;
    size_t nwords = buf[0];
    size_t nfills = buf[1];
    size  = buf[2];
    count = buf[3];
    words.resize(nwords);
    memcpy(words.data(),buf+4,nwords*sizeof(T));
    if (words.size() != nwords) {
        fprintf(stderr,"words.size() != nwords %zi %zi\n",words.size(),nwords);
        exit(1);
    }
    isFill.resize(nfills);
    memcpy(isFill.data(),buf + 4 + nwords, nfills*sizeof(T));
    if (isFill.size() != nfills) {
        fprintf(stderr,"isFill.size() != nfills %zi %zi\n",isFill.size(),nfills);
        exit(1);
    }

    firstActiveWord();
}

template <class T>
void BitVector<T>::firstActiveWord() {
    activeWordIdx   = 0;
    activeWordStart = 0;
    activeWordEnd   = nbits;
    activeWordType  = LITERAL;
    if (isFill[0] & (T)1) {
        if (words[0] >> (nbits-1))
            activeWordType = ONEFILL;
        else
            activeWordType = ZEROFILL;
        activeWordEnd = words[0] << shiftby;
    }
}

template <class T>
void BitVector<T>::nextActiveWord() {
    activeWordIdx++;
    activeWordStart = activeWordEnd;
    if (isFill[activeWordIdx >> shiftby] & ((T)1 << (activeWordIdx & (nbits-1)))) {
        if (words[activeWordIdx] >> (nbits-1))
            activeWordType = ONEFILL;
        else
            activeWordType = ZEROFILL;
        activeWordEnd += (words[activeWordIdx] << shiftby);
    }
    else {
        activeWordType  = LITERAL;
        activeWordEnd += nbits;
    }
}

// serialize the bitvector and return the number of words used
template <class T>
size_t BitVector<T>::dump(T **buf) {
    size_t nwords = 4 + words.size() + isFill.size();
    *buf = (T *) malloc(sizeof(T)*nwords);
    if (*buf == NULL) {
        fprintf(stderr,"failed to allocate %zi bytes\n",sizeof(T)*nwords);
        exit(4);
    }
    (*buf)[0] = (T) words.size();
    (*buf)[1] = (T) isFill.size();
    (*buf)[2] = (T) size;
    (*buf)[3] = (T) count;
    memcpy(*buf + 4,                 words.data(),  words.size()*sizeof(T));
    memcpy(*buf + 4 + words.size(), isFill.data(), isFill.size()*sizeof(T));
    return nwords;
}

// count the number of set bits in a word
template <class T>
int popcount(T v) {
    v = v - ((v >> 1) & (T)~(T)0/3);                           // temp
    v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);      // temp
    v = (v + (v >> 4)) & (T)~(T)0/255*15;                      // temp
    return (v * ((T)~(T)0/255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
}

// appends a word shaped uncompressed bitvector
template <class T>
void BitVector<T>::appendWord(T word) {
    if (word == (T)0) { // 0-fill
        if (activeWordType == ZEROFILL) // extends previous 0-fill
            words[activeWordIdx]++;
        else { // append a 0-fill
            int mod = words.size() & (nbits-1);
            if (mod)
                isFill[isFill.size()-1] |= (T)1 << mod;
            else
                isFill.push_back((T)1);
            words.push_back((T)1);
            activeWordType = ZEROFILL;
        }
    }
    else if (word == (T)~(T)0) { // 1-fill
        if (activeWordType == ONEFILL) // extends previous 1-fill
            words[activeWordIdx]++;
        else { // append a 1-fill
            int mod = words.size() & (nbits-1);
            if (mod)
                isFill[isFill.size()-1] |= (T)1 << mod;
            else
                isFill.push_back((T)1);
            words.push_back(onefill1); // ((T)1 << (nbits-1)) | (T)1
            activeWordType = ONEFILL;
            count += nbits;
        }
    }
    else { // literal
        if ((words.size() & (nbits-1)) == 0) // need another word in isFill
            isFill.push_back((T)0);
        words.push_back(word);
        activeWordType = LITERAL;
        count += popcount(word);
    }
    activeWordIdx   = words.size()-1;
    activeWordStart = size;
    size           += nbits;
    activeWordEnd   = size;
}

template <class T>
void BitVector<T>::seek(size_t wordStart) {
    if ((activeWordStart <= wordStart) && (wordStart < activeWordEnd)) return; // already here
    while (activeWordEnd <= wordStart) { // seek forward
        activeWordIdx++;
        activeWordStart = activeWordEnd;
        // check if this word is a fill
        if (isFill[activeWordIdx >> shiftby] & ((T)1 << (activeWordIdx & (nbits-1))))
            activeWordEnd += (words[activeWordIdx] << shiftby);
        else
            activeWordEnd += nbits;
    }
    while (activeWordStart > wordStart) { // seek backward
        activeWordIdx--;
        activeWordEnd = activeWordStart;
        // check if this word is a fill
        if (isFill[activeWordIdx >> shiftby] & ((T)1 << (activeWordIdx & (nbits-1))))
            activeWordStart -= (words[activeWordIdx] << shiftby);
        else
            activeWordStart -= nbits;
    }
    activeWordType = LITERAL;
    if (isFill[activeWordIdx >> shiftby] & ((T)1 << (activeWordIdx & (nbits-1)))) {
        if (words[activeWordIdx] >> (nbits-1))
            activeWordType = ONEFILL;
        else
            activeWordType = ZEROFILL;
    }
}

// fills a word shaped uncompressed bitvector starting at wordStart
template <class T>
void BitVector<T>::inflateWord(T *word, size_t wordStart) {
    seek(wordStart);
    if (activeWordType == LITERAL)
        *word = words[activeWordIdx];
    else if (activeWordType == ONEFILL)
        *word = (T)~(T)0;
    else
        *word = (T)0;
}

// extends the bitvector by length 0 bits (assumes previous bit was a 1)
template <class T>
void BitVector<T>::appendFill0(size_t length) {
    if (activeWordType == LITERAL) { // extend current LITERAL word
        T remainingBits = nbits - (size - activeWordStart);
        size += length;
        if (remainingBits >= length)
            return;
        // still more to do
        length -= remainingBits;
        if (length) activeWordStart += nbits;
    }
    else if (activeWordType == ONEFILL)
        size += length;
    else if (activeWordType == ZEROFILL) {
        // extend current zerofill
        size+=length;
        T nfills = (T)length >> shiftby;
        if (nfills) {
            words[words.size()-1] += nfills;
            length &= nbits-1;
        }
        // fprintf(stderr,"appendFill0 called when activeWordType == ZEROFILL!\n");
        // exit(5);
    }
    else // for the first word
        size += length;

    T nfills = (T)length >> shiftby;
    if (nfills) { // append a 0-fill
        // set a bit in isFill
        int mod = words.size() & (nbits-1);
        if (mod)
            isFill[isFill.size()-1] |= ((T)1 << mod);
        else
            isFill.push_back((T)1);
        words.push_back(nfills);
        activeWordIdx    = words.size()-1;
        activeWordType   = ZEROFILL;
        activeWordStart += (nfills << shiftby);
        length &= nbits-1;
    }
    if (length > 0) {
        if ((words.size() & (nbits-1)) == 0) isFill.push_back((T)0);
        words.push_back((T)0);
        activeWordType   = LITERAL;
        activeWordIdx    = words.size()-1;
    }
}
// extends the bitvector by length 1 bits (assumes previous bit was a 0)
template <class T>
void BitVector<T>::appendFill1(size_t length) {
    count += length;
    if (activeWordType == LITERAL) { // extend current LITERAL word
        int usedBits = size - activeWordStart;
        int remainingBits = nbits - usedBits;
        size += length;
        if (remainingBits) {
            if (length < remainingBits) {
                words[activeWordIdx] |= (((T)1 << length) - 1) << remainingBits - length;
                return;
            }
            else
                words[activeWordIdx] |= ((T)1 << remainingBits) - 1;
            length -= remainingBits;
        }
        if (length) activeWordStart += nbits;
    }
    else if (activeWordType == ZEROFILL) {
        size += length;
    }
    else if (activeWordType == ONEFILL) {
        // extend current onefill
        size+=length;
        T nfills = (T)length >> shiftby;
        if (nfills) {
            words[words.size()-1] += nfills;
            length &= nbits-1;
        }
        // fprintf(stderr,"appendFill1 called with activeWordType == ONEFILL\n");
        // exit(4);
    }
    else // for the first word
        size += length;

    T nfills = ((T)length >> shiftby);
    if (nfills) { // append a 1-fill
        // set a bit in isFill
        int mod = words.size() & (nbits-1);
        if (mod)
            isFill[isFill.size()-1] |= ((T)1 << mod);
        else
            isFill.push_back((T)1);
        words.push_back(nfills | ((T)1 << (nbits-1)));
        activeWordIdx    = words.size()-1;
        activeWordType   = ONEFILL;
        activeWordStart += (nfills << shiftby);
        length &= nbits-1;
    }
    if (length > 0) {
        // set length bits to 1 in literal word
        if ((words.size() & (nbits-1)) == 0) isFill.push_back((T)0);
        words.push_back((T)~(T)0 << (nbits - (int)length));
        activeWordType = LITERAL;
        activeWordIdx = words.size()-1;
    }
}

template <class T>
BitVector<T>* BitVector<T>::operator&(BitVector<T>* rhs) {
    BitVector<T> *res = new BitVector<T>();
    // first ensure that rhs is the same size.
    if (rhs->size != size) {
        fprintf(stderr,"rhs->size != size\n");
        return res;
    }
    firstActiveWord();
    rhs->firstActiveWord();
    while (activeWordIdx<words.size() and rhs->activeWordIdx < rhs->words.size()) {
        // advance until words overlap
        while (activeWordEnd <= rhs->activeWordStart) nextActiveWord();
        while (rhs->activeWordEnd <= activeWordStart) rhs->nextActiveWord();
        // compare overlapping words
        if (activeWordType == LITERAL) {
            if (rhs->activeWordType == LITERAL)
                res->appendWord(words[activeWordIdx] & rhs->words[rhs->activeWordIdx]);
            else if (rhs->activeWordType == ONEFILL)
                res->appendWord(words[activeWordIdx]);
            else
                res->appendFill0(rhs->activeWordEnd - activeWordStart);
            nextActiveWord();
        }
        else if (activeWordType == ONEFILL) {
            if (rhs->activeWordType == LITERAL) {
                res->appendWord(rhs->words[rhs->activeWordIdx]);
                rhs->nextActiveWord();
            }
            else if (rhs->activeWordType == ONEFILL) {
                if (activeWordEnd <= rhs->activeWordEnd) {
                    res->appendFill1(activeWordEnd - res->size);
                    nextActiveWord();
                }
                else {
                    res->appendFill1(rhs->activeWordEnd - res->size);
                    rhs->nextActiveWord();
                }
            }
            else {
                res->appendFill0(rhs->activeWordEnd - res->size);
                rhs->nextActiveWord();
            }
        }
        else { // ZEROFILL
            if (activeWordEnd <= rhs->activeWordEnd) {
                res->appendFill0(activeWordEnd - res->size);
                nextActiveWord();
            }
            else {
                res->appendFill0(rhs->activeWordEnd - res->size);
                rhs->nextActiveWord();
            }
        }
    }
    return res;
}

template <class T>
BitVector<T>* BitVector<T>::operator|(BitVector<T>* rhs) {
    BitVector<T> *res = new BitVector<T>();
    // first ensure that rhs is the same size.
    if (rhs->size != size) {
        fprintf(stderr,"rhs->size != size\n");
        return res;
    }
    firstActiveWord();
    rhs->firstActiveWord();
    while (activeWordIdx<words.size() and rhs->activeWordIdx < rhs->words.size()) {
        // advance until words overlap
        while (activeWordEnd <= rhs->activeWordStart) nextActiveWord();
        while (rhs->activeWordEnd <= activeWordStart) rhs->nextActiveWord();
        // compare overlapping words
        if (activeWordType == LITERAL) {
            if (rhs->activeWordType == LITERAL)
                res->appendWord(words[activeWordIdx] | rhs->words[rhs->activeWordIdx]);
            else if (rhs->activeWordType == ZEROFILL)
                res->appendWord(words[activeWordIdx]);
            else
                res->appendFill1(rhs->activeWordEnd - activeWordStart);
            nextActiveWord();
        }
        else if (activeWordType == ZEROFILL) {
            if (rhs->activeWordType == LITERAL) {
                res->appendWord(rhs->words[rhs->activeWordIdx]);
                rhs->nextActiveWord();
            }
            else if (rhs->activeWordType == ZEROFILL) {
                if (activeWordEnd <= rhs->activeWordEnd) {
                    res->appendFill0(activeWordEnd - res->size);
                    nextActiveWord();
                }
                else {
                    res->appendFill0(rhs->activeWordEnd - res->size);
                    rhs->nextActiveWord();
                }
            }
            else {
                res->appendFill1(rhs->activeWordEnd - res->size);
                rhs->nextActiveWord();
            }
        }
        else { // ONEFILL
            res->appendFill1(activeWordEnd - res->size);
            nextActiveWord();
        }
    }
    return res;
}

template <class T>
void BitVector<T>::operator&=(BitVector<T>* rhs) {
    BitVector<T>* res = *this & rhs;
    words.swap(res->words);
    isFill.swap(res->isFill);
    count = res->count;
}

template <class T>
void BitVector<T>::operator|=(BitVector<T>* rhs) {
    BitVector<T>* res = *this | rhs;
    words.swap(res->words);
    isFill.swap(res->isFill);
    count = res->count;
}

template <class T>
void BitVector<T>::flip() {
    // iterate over the words
    // flip fill words by toggling MSB
    // flip literal words with XOR ~0
    T allones = ~(T)0;
    T one2zero = allones >> 1;
    T zero2one = ~one2zero;
    for(size_t i=0;i<words.size();i++)
        if(isFill[i >> shiftby] & ((T)1 << (i & (nbits-1))))
            if (words[i] >> (nbits-1)) // one-fill
                words[i] &= one2zero;
            else
                words[i] |= zero2one;
        else
            words[i] ^= allones;
}

#endif
