#include "bvec32.h"
#include <algorithm>

#define DEBUG false

void print_binary(vector<uint32_t>& v) {
    for (int i = 0; i < v.size(); i++) {
        for (int b = 31; b >= 0; b--)
            printf("%d", v[i] & (1 << b) ? 1 : 0);
		printf(" ");
    }
    printf("\n");
}

void xor_vectors(const vector<uint32_t>& a, const vector<uint32_t>& b,
    vector<uint32_t>& c) {
    for (int i = 0; i < a.size(); i++) {
        c.push_back(a[i] ^ b[i]);
    }
}

void debug_binary(const char* head, vector<uint32_t>& v) {
    if (!DEBUG) return;
    printf("%-20s ", head);
    print_binary(v);
}

using namespace std;
int main(int argc, char *argv[]) {
	// simulate random bitvectors with up to 1 billion points
	// 16mers - up to 2^32 - 1
	// 32mers - up to 2^64 - 1
	srand(time(NULL));
	uint32_t n = 1000;
	if (argc > 1) 
		n = atoi(argv[1]);

	vector<uint32_t> rand32;
	for(int i=0;i<n;i++) {
		rand32.push_back(rand() | ((rand() % 4) << 30));
	}
	sort(rand32.begin(),rand32.end());
	printf("max value selected: %u\n",rand32.back());
	vector<uint32_t>::iterator ii32 = rand32.begin();
	vector<uint32_t> rand32uniq;
	rand32uniq.push_back(*ii32);
	++ii32;
	while (ii32 != rand32.end()) {
		if (*ii32 != rand32uniq.back())
			rand32uniq.push_back(*ii32);
		++ii32;
	}
	printf("generated %zi distinct random uint32_t\n",rand32uniq.size());
	bvec *bv32 = new bvec(rand32uniq);
//    bv32->print();
	printf("cnt(): %u, bv32->bytes(): %u\n",bv32->cnt(),bv32->bytes());
	bv32->compress();
    vector<uint32_t> compressed = bv32->get_words();
//    bv32->print();
	printf("cnt(): %u, bv32->bytes(): %u\n",bv32->cnt(),bv32->bytes());
	bv32->decompress();
    vector<uint32_t> uncompressed = bv32->get_words();
//    bv32->print();
	printf("cnt(): %u, bv32->bytes(): %u\n",bv32->cnt(),bv32->bytes());
	delete bv32;
    cout << "Are they equal? " << (rand32uniq == uncompressed ? "yes" : "no") << endl;
    if (rand32uniq != uncompressed) {
        vector<uint32_t> v;
        xor_vectors(rand32uniq, uncompressed, v);
        debug_binary("Original", rand32uniq);
		debug_binary("Compressed", compressed);
        debug_binary("Uncompressed", uncompressed);
        debug_binary("XOR", v);
    }

	// vector<uint64_t> rand64;
	// for(int i=0;i<n;i++) {
	// 	rand64.push_back((uint64_t)rand() | ((uint64_t)rand() << 30) | ((uint64_t)rand() << 60));
	// }
	// sort(rand64.begin(),rand64.end());
	// printf("max value selected: %llu\n",rand64.back());
	// vector<uint64_t>::iterator ii64 = rand64.begin();
	// vector<uint64_t> rand64uniq;
	// rand64uniq.push_back(*ii64);
	// ++ii64;
	// while (ii64 != rand64.end()) {
	// 	if (*ii64 != rand64uniq.back())
	// 		rand64uniq.push_back(*ii64);
	// 	++ii64;
	// }
	// printf("generated %zi distinct random uint64_t\n",rand64uniq.size());
	// bvec *bv64 = new bvec(rand64uniq);
	// printf("bv64()->bytes() = %u\n",bv64->bytes());
	// delete bv64;
}