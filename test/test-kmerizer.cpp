#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Kmerizer Unit Tests"

#include <boost/test/unit_test.hpp>
#include "test.h"
#include "kmerizer/kmerizer.h"

Kmerizer * initKmerizer() {
    return new Kmerizer(1, 1, "/tmp", CANONICAL);
}

BOOST_AUTO_TEST_SUITE(Kmerizer);
    BOOST_AUTO_TEST_CASE(Instantiation) {
        BOOST_REQUIRE(initKmerizer() != NULL);
        cout << "1..1" << endl;
        cout << "ok 1" << endl;
    }
BOOST_AUTO_TEST_SUITE_END();
