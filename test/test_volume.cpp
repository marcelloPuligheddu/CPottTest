#include <gtest/gtest.h>
#include <random>
#include "CellularPotts.h"

TEST(CPM, VolumeConservation) {
    std::mt19937 gen(1234);
    CellularPotts cpm(32, 10, gen);

    for (int i = 0; i < 1000; ++i)
        cpm.monte_carlo_step();

    EXPECT_EQ(cpm.total_volume(), 32*32);
}
