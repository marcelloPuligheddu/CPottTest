#include <gtest/gtest.h>
#include <random>
#include "CellularPotts.h"

TEST(CPM, ZeroTemperatureRejectsPositiveDelta)
{
    std::mt19937 gen(1234);
    CellularPotts<16, 4> cpm (gen);
    cpm.set_temperature(0.0);

    EXPECT_FALSE(cpm.accept(5.0));
    EXPECT_TRUE(cpm.accept(-1.0));
}
