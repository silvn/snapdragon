#include "test.h"
#include "kmerizer/kmerizer.h"

int run_test_cases() {
    Kmerizer * km = new Kmerizer(1, 1, "/tmp", CANONICAL);
    return 0;
}

int main(int argc, char *argv[]) {
    return run_test_cases();
}