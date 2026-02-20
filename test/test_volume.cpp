#include <gtest/gtest.h>
#include <random>
#include "CellularPotts.h"

TEST(CPM, ApplyMoveUpdatesVolume) {
    std::mt19937 gen(1234);
    CellularPotts cpm<16, 4> (gen);

    int i = 5;
    int j = 6;    
    
    int before = cpm.total_volume();
    cpm.apply_move(i, j);
    EXPECT_EQ(cpm.total_volume(), before);
    
}
