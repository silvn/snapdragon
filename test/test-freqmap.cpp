#include "test.h"
#include "util/freqmap.h"

namespace {
    
class FrequencyMapTest : public ::testing::Test {
};

TEST(FrequencyMapTest, DoesInstantiation) {
    FrequencyMap *map = new FrequencyMap();
    EXPECT_TRUE(map != NULL);
}

} /* namespace */

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
