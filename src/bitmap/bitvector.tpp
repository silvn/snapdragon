void printBits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    for (i=size-1;i>=0;i--)
    {
        for (j=7;j>=0;j--)
        {
            byte = b[i] & (1<<j);
            byte >>= j;
           fprintf(stderr,"%u", byte);
        }
    }
}


template <class T>
BitVector<T>::BitVector() {
    count=0;
    size=0;
    // isFill.push_back(0);
    activeWordIdx = 0;
    activeWordStart = 0;
    activeWordEnd = 0;
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
    nbits = sizeof(T)*8;
    if (nbits == 64) shiftby = 6;
    else if (nbits == 32) shiftby = 5;
    else if (nbits == 16) shiftby = 4;
    else shiftby = 3;
    size_t nwords = buf[0];
    size_t nfills = buf[1];
    size = buf[2];
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
    activeWordIdx = 0;
    activeWordStart=0;
    activeWordEnd = (size < nbits) ? size-1 : nbits-1;
    activeWordType = LITERAL;
    if (isFill[0] & 1) {
        if (words[0] >> (nbits-1))
            activeWordType = ONEFILL;
        else
            activeWordType = ZEROFILL;
        activeWordEnd = (words[0] << shiftby) - 1;
    }
}

// DIY serialization
template <class T>
size_t BitVector<T>::dump(T **buf) {
    // allocate space in buf
    size_t dbytes = sizeof(T)*(4 + words.size() + isFill.size());
    *buf = (T *) malloc(dbytes);
    if (*buf == NULL) {
        fprintf(stderr,"failed to allocate %zi bytes\n",dbytes);
        return 0;
    }
    (*buf)[0] = (T) words.size();
    (*buf)[1] = (T) isFill.size();
    (*buf)[2] = (T) size;
    (*buf)[3] = (T) count;
    memcpy(*buf + 4, words.data(), words.size()*sizeof(T));
    memcpy(*buf + 4 + words.size(), isFill.data(), isFill.size()*sizeof(T));
    return 4 + words.size() + isFill.size();
}

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
    // fprintf(stderr,"appendWord() size %i activeWordIdx %i activeWordStart %i isFill.size() %i words.size() %i\n",size,activeWordIdx,activeWordStart,isFill.size(),words.size());
    activeWordEnd+=nbits;
    if (word == (T)0) { // 0-fill
        if (activeWordType == ZEROFILL)
            words[activeWordIdx]++;
        else {
            /* set calculated bit */;
            int mod = words.size() & (nbits-1);
            if (mod)
                isFill[isFill.size()-1] = isFill.back() | ((T)1 << mod);
            else
                isFill.push_back((T)1);
            words.push_back((T)1);
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
            int mod = words.size() & (nbits-1);
            if (mod)
                isFill[isFill.size()-1] = isFill.back() | ((T)1 << mod);
            else
                isFill.push_back((T)1);
            words.push_back(((T)1 << (nbits-1)) | (T)1);
            activeWordIdx = words.size()-1;
            activeWordType = ONEFILL;
            activeWordStart = size;
            count += nbits;
        }
    }
    else { // literal
        if ((words.size() & (nbits-1)) == 0)
            isFill.push_back((T)0);
        words.push_back(word);
        activeWordIdx = words.size()-1;
        activeWordType = LITERAL;
        activeWordStart = size;
        count += popcount(word);
    }
    size += nbits;
}

template <class T>
void BitVector<T>::seek(size_t wordStart) {
    if (activeWordStart == wordStart) return;
    while (activeWordEnd < wordStart) { // seek forward
        activeWordIdx++;
        activeWordStart = activeWordEnd+1;
        // check if this word is a fill
        if (isFill[activeWordIdx >> shiftby] & ((T)1 << (activeWordIdx & (nbits-1))))
            activeWordEnd += words[activeWordIdx] << shiftby;
        else
            activeWordEnd += nbits;
    }
    while (activeWordStart > wordStart) { // seek backward
        activeWordIdx--;
        activeWordEnd = activeWordStart-1;
        // check if this word is a fill
        if (isFill[activeWordIdx >> shiftby] & ((T)1 << (activeWordIdx & (nbits-1))))
            activeWordStart -= words[activeWordIdx] << shiftby;
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

// decompress nbits of the bitvector starting at wordStart
// for each set bit, copy v into A[bitpos]
template <class T>
void BitVector<T>::fillSetBits(size_t wordStart, uint32_t *A, uint32_t v) {
    seek(wordStart);
    if (activeWordType == LITERAL) {
        T bits = words[activeWordIdx];
        while (bits) {
            A[ffs(bits)-1] = v;
            bits &= bits-1;
        }
    }
    else if (activeWordType == ONEFILL)
        for(int i=0;i<nbits;i++)
            A[i]=v;
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
        activeWordStart += nbits;
    }
    else if (activeWordType == ONEFILL) {
        activeWordStart += words[activeWordIdx] << shiftby;
        size += length;
    }
    else if (activeWordType == ZEROFILL) {
        fprintf(stderr,"appendFill0 called when activeWordType == ZEROFILL!\n");
        exit(5);
    }
    else // for the first word
        size += length;

    T nfills = length >> shiftby;
    if (nfills) { // append a 0-fill
        // set a bit in isFill
        T mod = words.size() & (nbits-1);
        if (mod)
            isFill[isFill.size()-1] = isFill.back() | ((T)1 << mod);
        else
            isFill.push_back((T)1);
        words.push_back(nfills);
        activeWordIdx = words.size()-1;
        activeWordType = ZEROFILL;
        length &= nbits-1;
        activeWordStart += nbits * nfills;
    }
    if (length > 0) {
        words.push_back((T)0);
        if ((words.size() >> shiftby) == isFill.size()) isFill.push_back((T)0);
        activeWordType = LITERAL;
        activeWordIdx = words.size()-1;
    }
}
// extends the bitvector by length 1 bits (assumes previous bit was a 0)
template <class T>
void BitVector<T>::appendFill1(size_t length) {
    if (activeWordType == LITERAL) { // extend current LITERAL word
        T usedBits = size - activeWordStart;
        size += length;
        words[activeWordIdx] |= ((T)~(T)0) << usedBits; // fill the word with 1's
        T remainingBits = nbits - usedBits;
        if (remainingBits > length) { // flip some bits back to 0
            words[activeWordIdx] &= ((T)~(T)0) >> (remainingBits - length);
            return;
        }
        length -= remainingBits;
        activeWordStart += nbits;
    }
    else if (activeWordType == ZEROFILL) {
        activeWordStart += words[activeWordIdx] << shiftby;  
        size += length;
    }
    else if (activeWordType == ONEFILL) {
        fprintf(stderr,"appendFill1 called with activeWordType == ONEFILL\n");
        exit(4);
    }
    else // for the first word
        size += length;

    T nfills = length >> shiftby;
    if (nfills) { // append a 1-fill
        // set a bit in isFill
        T mod = words.size() & (nbits-1);
        if (mod)
            isFill[isFill.size()-1] = isFill.back() | ((T)1 << mod);
        else
            isFill.push_back((T)1);
        words.push_back(((T)1 << (nbits-1)) | nfills);
        activeWordIdx = words.size()-1;
        activeWordType = ONEFILL;
        length &= nbits-1;
        activeWordStart += nfills*nbits;
    }
    if (length > 0) {
        // set length bits to 1 in literal word
        words.push_back((T)~(T)0 >> (nbits - length));
        if ((words.size() >> shiftby) == isFill.size()) isFill.push_back((T)0);
        activeWordType = LITERAL;
        activeWordIdx = words.size()-1;
    }
}
