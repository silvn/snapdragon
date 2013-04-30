#ifndef SNAPDRAGON_BVEC_RIWAH_H
#define SNAPDRAGON_BVEC_RIWAH_H

#include <stdint.h>

#define DEBUG false
#define LITERAL_SIZE 32

typedef uint32_t word_t;

/*
    Recursively Indexed Word Aligned Hybrid (RIWAH) compression
    A bitvector may be compressed into a mix of 32 bit literal and fill words.
    Instead of coding the word type within the word in the first 2 bits, as in Wu et al's WAH
    we use external data structures to track the types of words.
    
    word types:
        literal: all 32 bits are uncompressed
        0-fill: value holds the length of a contiguous run of 0's (multiple of 32)
        1-fill: value holds the length of a contiguous run of 1's (multiple of 32)

        dirty fills:
            A fill word can sometimes be merged with the next literal word (if space allows).
            The most significant bits hold the position(s) of the differing bit(s) in the next 32.
            The remaining bits in the word hold the length of the fill.

            The 0-fill and 1-fill word types can be thought of as dirty fills with 0 different
            bits in the next word. So there are a total of 14 types of fill words
            [01]-dfill-0:  0 pos bits, 32 length bits - maxfill 2^(32+5) = 137438953472
            [01]-dfill-1:  5 pos bits, 27 length bits - maxfill 2^(27+5) =   4294967296
            [01]-dfill-2: 10 pos bits, 22 length bits - maxfill 2^(22+5) =    134217728
            [01]-dfill-3: 15 pos bits, 17 length bits - maxfill 2^(17+5) =      4194304
            [01]-dfill-4: 20 pos bits, 12 length bits - maxfill 2^(12+5) =       131072
            [01]-dfill-5: 25 pos bits,  7 length bits - maxfill 2^(7+5)  =         4096
            [01]-dfill-6: 30 pos bits,  2 length bits - maxfill 2^(2+5)  =          128

        merged lonely bit literals:
            a sequence of adjacent literal words may be compressed using the same concepts.
            If 5 bits can indicate the position of a bit, we could compress up to 6 adjacent lonely bit literals
            into one 32 bit word with 2 bits left over. I think 6 adjacent lonely bit literals is unlikely, so lets do 5.
            The first 2 bits indicate the number of literals-2.
            The next 5 bits each indicate whether a literal needs to be flipped.
            The remaining 25 bits hold the positions of the lonely bit in each literal.

    We need to have an index for each of the 15 types of compressed words in the bit vector. The index may be just an array
    of offsets or another RIWAH compressed bit vector.  We can let the new RIWAH decide for us. Some formula will be used
    to determine whether to use an array or a compressed bit vector.
    
    considerations:
    bitwise logical operations
        keep track of the next compressed word and its type.
        need to decompress and compress words in situ
        The word indexes need to be updated on the fly as well ( tricky for &=, |=, and ^=. create new ones and swap?)

    searching/iterating - use the active word trick to avoid having to do a linear scan
        find(x)
        nextOne(x)

    building/modifying
        constructors
        appendFill()
        setBit() - might be expensive - can we use the word type indexes to locate the word in < O(n)?

*/
#endif // SNAPDRAGON_BVEC_RIWAH_H