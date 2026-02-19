double E0 = cpm.total_energy();
double dH = cpm.attempt_copy_once_and_return_dH();
double E1 = cpm.total_energy();

EXPECT_NEAR(E1 - E0, dH, 1e-12);
