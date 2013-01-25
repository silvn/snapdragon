#include "bvec.h"

int main(int argc, char *argv[]) {
	vector<uint64_t> vec1,vec2;
	vec1.push_back(0);
	vec1.push_back(100);
	vec2.push_back(0);
	vec2.push_back(10);
	bvec *bv1 = new bvec(vec1);
	bvec *bv2 = new bvec(vec2);
	printf("bv1->cnt()=%llu\n",bv1->cnt());
	printf("bv2->cnt()=%llu\n",bv2->cnt());
	bvec *bvu = *bv1 | *bv2;
	printf("bvu->cnt()=%llu\n",bvu->cnt());
	bvec *bvi = *bv1 & *bv2;
	printf("bvi->cnt()=%llu\n",bvi->cnt());
	return 0;
}