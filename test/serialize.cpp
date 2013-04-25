#include "bvec/bvec.h"
#include "test.h"
#include <algorithm>

BitVector* random_bvec(int n) {
    // simulate random bitvectors with up to 1 billion points
    // 16mers - up to 2^32 - 1
    // 32mers - up to 2^64 - 1
    srand(time(NULL));

    vector<uint32_t> rand32;
    for(int i=0;i<n;i++) {
        rand32.push_back(rand() | ((rand() % 4) << 30));
    }
    sort(rand32.begin(),rand32.end());
    printf("max value selected: %u\n",rand32.back());
    vector<uint32_t>::iterator ii32 = rand32.begin();
    vector<uint32_t> uniq;
    uniq.push_back(*ii32);
    ++ii32;
    while (ii32 != rand32.end()) {
        if (*ii32 != uniq.back())
            uniq.push_back(*ii32);
        ++ii32;
    }
    printf("generated %zi distinct random uint32_t\n",uniq.size());
    return new BitVector(uniq);
}

using namespace std;
int main(int argc, char *argv[]) {
    short n = 1000;
    if (argc > 1) n = atoi(argv[1]);

    BitVector * original = random_bvec(n);
    BitVector * deserialized = new BitVector();

    char filename[200];
    std::tmpnam(filename);
    cout << "Writing to file " << filename << endl;

    BitVector::save(*original, filename);
    BitVector::restore(*deserialized, filename);

    bool match = original->equals(*deserialized); // original == deserialized
    if (!match) {
        original->print();
        deserialized->print();
        vector<uint32_t> orig = original->getWords();
        vector<uint32_t> deser = deserialized->getWords();
        debug_binary("Original", orig);
        debug_binary("Deserialized", deser);
        vector<uint32_t> v;
        xor_vectors(orig, deser, v);
        debug_binary("XOR", v);
        return 1;
    }
    return 0;
}