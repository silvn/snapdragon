#include "bvec.h"

int main(int argc, char *argv[]) {
	vector<uint64_t> vec64;
	vec64.push_back(0);
	vec64.push_back(100);
	bvec *bv = new bvec(vec64);
	printf("bv->cnt()=%llu\n",bv->cnt());
	return 0;
}