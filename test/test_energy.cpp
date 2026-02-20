#include <gtest/gtest.h>
#include <random>
#include "CellularPotts.h"

TEST(CPM, DeltaMatchesFullEnergyDifference)
{
    std::mt19937 gen(1234);
    CellularPotts<16, 4> cpm (gen);

    int i = 5;
    int j = 6;

    double E0 = cpm.total_energy();
    double dH = cpm.compute_delta(i, j);

    cpm.apply_move(i, j);
    double E1 = cpm.total_energy();

    EXPECT_NEAR(E1 - E0, dH, 1e-12);
}
