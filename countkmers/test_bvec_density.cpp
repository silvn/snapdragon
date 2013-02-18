#include "bvec32.h"
#include <algorithm>

#define TERM_RESET "\e[m"
#define TERM_BOLD  "\e[1m"
#define TERM_RED   "\e[31m"
#define TERM_GREEN "\e[32m"

#define TEST_VEC_LENGTH 4
const uint32_t TEST_CASES[][TEST_VEC_LENGTH] = {
    { 161143029,  1089249216, 1093915745, 1157787156 },
    { 3085261341, 3247604874, 3456053010, 3773962658 },
    { 13431029,   377786868,  1463271703, 1501810865 },
    { 1819233732, 2235579467, 2625766209, 4061912144 }
};

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
    printf("%-20s ", head);
    print_binary(v);
}

void run_test_cases() {
    int num_cases = sizeof(TEST_CASES) / (TEST_VEC_LENGTH * sizeof(uint32_t));
    for (int i = 0; i < num_cases; i++) {
        vector<uint32_t> v;
        for (int j = 0; j < TEST_VEC_LENGTH; j++)
            v.push_back(TEST_CASES[i][j]);
        bvec32 *bv = new bvec32(v);
        bv->compress();
        bv->decompress();
        vector<uint32_t> uncomp = bv->get_words();
        bool result = v == uncomp;
        printf(""TERM_BOLD"Test %d"TERM_RESET": %s\n", i, result
            ? ""TERM_GREEN"Pass"TERM_RESET"" : ""TERM_RED"Fail"TERM_RESET"");
        if (!result) {
            vector<uint32_t> va;
            printf("%-20s ", "Original");
            for (int j = 0; j < v.size(); j++) printf("%32u ", v[j]);
            cout << endl;
            xor_vectors(v, uncomp, va);
            debug_binary("Original (bin)", v);
            debug_binary("Uncompressed", uncomp);
            debug_binary("XOR", va);
        }
    }
}

using namespace std;
int main(int argc, char *argv[]) {
    run_test_cases();
    return 0;
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
	bvec32 *bv32 = new bvec32(rand32uniq);
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
	// bvec32 *bv64 = new bvec32(rand64uniq);
	// printf("bv64()->bytes() = %u\n",bv64->bytes());
	// delete bv64;
}