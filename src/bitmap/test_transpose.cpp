#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <bitset>
// #include "bitmap.h"
using namespace std;

void print_binary(unsigned long A[64]) {
    for (int i = 0; i < 64; i++) {
        cout << bitset<64>(A[i]).to_string() << " " << i << endl;
    }
}
void transpose64(unsigned long A[64]) {
    int j, k;
    unsigned long m, t;
    
    m = 0x00000000FFFFFFFF;
    for (j = 32; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 64; k = (k + j + 1) & ~j) {
            t = (A[k] ^ (A[k+j] >> j)) & m;
            A[k] = A[k] ^ t;
            A[k+j] = A[k+j] ^ (t << j);
        }
    }
}

int main() {
    unsigned long A[64],B[64];
    for(unsigned long i=0;i<64;i++) {
        A[i] = 0;//xFFFFFFFF00000000;
        for(int j=0;j<8;j++) {
            A[i] |= i << 8*j; //((unsigned long)(rand() & 255)) << (8*j);
        }
    }
    print_binary(A);
    for(int i =0; i < 10000001; i++) {
        transpose64(A);
    }
    print_binary(A);
    return 0;
}