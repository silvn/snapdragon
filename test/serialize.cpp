#include "bvec32.h"
#include "test.h"
#include <algorithm>

bvec32* random_bvec(int n) {
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
    return new bvec32(uniq);
}

using namespace std;
int main(int argc, char *argv[]) {
	short n = 1000;
	if (argc > 1) 
		n = atoi(argv[1]);

	bvec32 * original = random_bvec(n);
    bvec32 * deserialized = new bvec32();
    
    char filename[200];
    std::tmpnam(filename);
    cout << "Writing to file " << filename << endl;
    
    // save_to_file(*original, filename);
    // restore_from_file(*deserialized, filename);
    
    print_test("serialize", original == deserialized);
    return 0;
}