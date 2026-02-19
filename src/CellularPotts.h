#include <vector>
#include <random>
#include <cmath>
#include <iostream>

class CellularPotts {
public:
    CellularPotts(int L, int num_cells, std::mt19937 gen_);
    void monte_carlo_step();
    double total_energy() const;
    int total_volume() const;
private:
    int L;
    std::vector<int> lattice;
    std::vector<int> volume;
    std::vector<int> V0;
    std::vector<std::vector<double>> J;
    double lambdaV;
    double temperature;

    std::mt19937 gen;
    std::uniform_int_distribution<int> dist_site;
    std::uniform_real_distribution<double> dist_prob;

    inline int idx(int x, int y) const ;
    int neighbor(int x, int y) ;
    double adhesion_energy(int i, int new_sigma) ;
    double volume_energy(int sigma_old, int sigma_new) ;
    double attempt_copy() ;
};
