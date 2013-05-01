#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "BitVector Serialization Unit Tests"

#include <boost/test/unit_test.hpp>
#include "bvec/bvec.h"
#include "test.h"
#include <algorithm>

#define TEST_VEC_LENGTH 4
const uint32_t TEST_CASES[][TEST_VEC_LENGTH] = {
    { 161143029,  1089249216, 1093915745, 1157787156 },
    { 3085261341, 3247604874, 3456053010, 3773962658 },
    { 13431029,   377786868,  1463271703, 1501810865 },
    { 1819233732, 2235579467, 2625766209, 4061912144 }
};

using namespace std;
BOOST_AUTO_TEST_SUITE(BitVectorDensity);
    BOOST_AUTO_TEST_CASE(RoundTrip) {
        int num_cases = sizeof(TEST_CASES) /
            (TEST_VEC_LENGTH * sizeof(uint32_t));
        int retValue = 0;
        for (int i = 0; i < num_cases; i++) {
            vector<uint32_t> v;
            for (int j = 0; j < TEST_VEC_LENGTH; j++)
                v.push_back(TEST_CASES[i][j]);
            BitVector *bv = new BitVector(v);
            bv->compress();
            bv->decompress();
            vector<uint32_t> uncomp = bv->getWords();
            BOOST_CHECK(v == uncomp);
        }
    }
}
