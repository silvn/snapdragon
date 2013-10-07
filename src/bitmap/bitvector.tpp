template <class T>
BitVector<T>::BitVector() {
    count=0;
    size=0;
    // words.reserve(30000); // this is a hack because the activeWord is a regular pointer and gets messed up if words is resized
    // isFill.reserve(32);
    isFill.push_back(0);
    activeWordIdx = 0;
    activeWordStart = 0;
    activeWordType = NONEOFTHEABOVE;
    nbits = sizeof(T)*8;
    if (nbits == 64) shiftby = 6;
    else if (nbits == 32) shiftby = 5;
    else if (nbits == 16) shiftby = 4;
    else shiftby = 3;
}

// constructor - given a previously dumped BitVector
template <class T>
BitVector<T>::BitVector(T *buf) {
    T nwords = buf[0];
    T nfills = buf[1];
    size = buf[2];
    count = buf[3];
    words.reserve(nwords);
    memcpy(words.data(),buf+4,nwords*sizeof(T));
    isFill.reserve(nfills);
    memcpy(isFill.data(),buf + 4 + nwords, nfills*sizeof(T));
    activeWordIdx = 0;
    activeWordStart=0;
    activeWordType = LITERAL;
    if (isFill[0] & 1) {
        if (words[0] & 1 << (nbits-1))
            activeWordType = ONEFILL;
        else
            activeWordType = ZEROFILL;
    }
    nbits = sizeof(T)*8;
    if (nbits == 64) shiftby = 6;
    else if (nbits == 32) shiftby = 5;
    else if (nbits == 16) shiftby = 4;
    else shiftby = 3;
}

// DIY serialization
template <class T>
size_t BitVector<T>::dump(T **buf) {
    // allocate space in buf
    size_t dbytes = sizeof(T)*(4 + words.size() + isFill.size());
    *buf = malloc(dbytes);
    if (*buf == NULL) {
        fprintf(stderr,"failed to allocate %zi bytes\n",dbytes);
        return 0;
    }
    (*buf)[0] = (T)words.size();
    (*buf)[1] = (T)isFill.size();
    (*buf)[2] = size;
    (*buf)[3] = count;
    memcpy(*buf + 4, words.data(), words.size()*sizeof(T));
    memcpy(*buf + 4 + words.size(), isFill.data(), isFill.size()*sizeof(T));
    return dbytes;
}

template <class T>
T popcount(T v) {
    v = v - ((v >> 1) & (T)~(T)0/3);                           // temp
    v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);      // temp
    v = (v + (v >> 4)) & (T)~(T)0/255*15;                      // temp
    return (T)(v * ((T)~(T)0/255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
}

// appends a word shaped uncompressed bitvector
template <class T>
void BitVector<T>::appendWord(T word) {
    size += nbits;
    if (word == (T)0) { // 0-fill
        if (activeWordType == ZEROFILL)
            words[activeWordIdx]++;
        else {
            /* set calculated bit */;
            T mod = words.size() & (nbits-1);
            if (mod)
                isFill[isFill.size()-1] = isFill.back() | (1 << (mod-1));
            else
                isFill.push_back(1);
            words.push_back(1);
            activeWordIdx = words.size()-1;
            activeWordType = ZEROFILL;
            activeWordStart = size;
        }
    }
    else if (word == (T)~(T)0) { // 1-fill
        if (activeWordType == ONEFILL)
            words[activeWordIdx]++;
        else {
            /* set calculated bit */
            T mod = words.size() & (nbits-1);
            if (mod)
                isFill[isFill.size()-1] = isFill.back() | (1 << (mod-1));
            else
                isFill.push_back(1);
            words.push_back((1 << (nbits-1)) | 1);
            activeWordIdx = words.size()-1;
            activeWordType = ONEFILL;
            activeWordStart = size;
            count += nbits;
        }
    }
    else { // literal
        if (words.size() & (nbits-1) == 0)
            isFill.push_back(0);
        words.push_back(word);
        activeWordIdx = words.size()-1;
        activeWordType = LITERAL;
        activeWordStart = size;
        count += popcount(word);
    }
}

// fills a word shaped uncompressed bitvector starting at bit_pos
// returns false if end of vector
template <class T>
bool BitVector<T>::nextWord(T *word, T *bit_pos) {
    if (activeWordIdx == words.size()-1) return false;

    *bit_pos += nbits;

    if (activeWordType == LITERAL)
        *word = words[activeWordIdx];
    else {
        if (activeWordType == ZEROFILL)
            *word = (T)0;
        else
            *word = (T)~(T)0;
        // check if we are not yet at the end of the word
        if (*bit_pos < activeWordStart + words[activeWordIdx] & (~(T)0 >> 1))
            return true;
    }
    activeWordIdx++;
    activeWordType = LITERAL;
    activeWordStart = *bit_pos;
    // is the current word a fill word?
    T offset = (activeWordIdx - words.start()) >> shiftby;
    T bit = 1 << ((activeWordIdx - words.start()) & (nbits-1));
    if (isFill[offset] & bit)
        activeWordType = (words[activeWordIdx] & (1 << (nbits-1))) ? ONEFILL : ZEROFILL;
    return true;
}

// extends the bitvector by length 0 bits (assumes previous bit was a 1)
template <class T>
void BitVector<T>::appendFill0(T length) {
    // fprintf(stderr,"appendFill0(%llu) activeWordType: %i, size: %llu, activeWordStart: %llu\n",length,activeWordType,size,activeWordStart);
    if (activeWordType == LITERAL) { // extend current LITERAL word
        T remainingBits = nbits - (size - activeWordStart);
        size += length;
        if (remainingBits >= length)
            return;
        length -= remainingBits;
        activeWordStart += nbits;
    }
    else if (activeWordType == ONEFILL) {
        activeWordStart += nbits * (words[activeWordIdx] & ((T)~(T)0 >> 1));  
        size += length;
    }
    
    T nfills = length >> shiftby;
    // fprintf(stderr,"zzzlength: %llu, nfills: %llu\n",length,nfills);
    if (nfills) { // append a 0-fill
        // fprintf(stderr,"before\n");
        words.push_back(nfills);
        // fprintf(stderr,"after\n");
        activeWordIdx = words.size()-1;
        activeWordType = ZEROFILL;
        // set a bit in isFill
        T mod = words.size() & (nbits-1);
        // fprintf(stderr,"mod: %llu\n",mod);
        if (mod)
            isFill[isFill.size()-1] = isFill.back() | ((T)1 << (mod-1));
        else
            isFill.push_back((T)1);
        length &= nbits-1;
        if (length == 0) return;
        activeWordStart += nbits * nfills;
    }
    // fprintf(stderr,"length: %llu\n",length);
    if (length > 0) {
        words.push_back((T)0);
        if ((words.size() >> shiftby) == isFill.size()) isFill.push_back((T)0);
        activeWordType = LITERAL;
        activeWordIdx = words.size()-1;
    }
}
// extends the bitvector by length 1 bits (assumes previous bit was a 0)
template <class T>
void BitVector<T>::appendFill1(T length) {
    if (activeWordType == LITERAL) { // extend current LITERAL word
        T usedBits = size - activeWordStart;
        size += length;
        words[activeWordIdx] |= ((T)~(T)0) >> usedBits; // fill the word with 1's
        T remainingBits = nbits - usedBits;
        if (remainingBits > length) { // flip some bits back to 0
            words[activeWordIdx] &= ((T)~(T)0) << (remainingBits - length);
            return;
        }
        length -= remainingBits;
        activeWordStart += nbits;
    }
    else if (activeWordType == ZEROFILL) {
        activeWordStart += words[activeWordIdx] * nbits;  
        size += length;
    }
    T nfills = length >> shiftby;
    if (nfills) { // append a 1-fill
        words.push_back(((T)1 << (nbits-1)) | nfills);
        activeWordIdx = words.size()-1;
        activeWordType = ONEFILL;
        // set a bit in isFill
        T mod = words.size() & (nbits-1);
        if (mod)
            isFill[isFill.size()-1] = isFill.back() | ((T)1 << (mod-1));
        else
            isFill.push_back((T)1);
        length &= nbits-1;
        if (length == 0) return;
        activeWordStart += nfills*nbits;
    }
    if (length > 0) {
        // set length bits to 1 in literal word
        words.push_back((T)~(T)0 << (nbits - length));
        if ((words.size() >> shiftby) == isFill.size()) isFill.push_back((T)0);
        activeWordType = LITERAL;
        activeWordIdx = words.size()-1;
    }
}
