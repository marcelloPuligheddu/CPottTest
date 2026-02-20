#include <vector>
#include <random>
#include <cmath>
#include <iostream>

template< int L, int K >
class CellularPotts {
public:
    CellularPotts(int L, int num_cells, std::mt19937& gen);

    void monte_carlo_step();
    double compute_delta(int i, int j) const;
    void apply_move(int i, int j);
    bool accept(double dH);
    bool attempt_copy();
    double total_energy() const;
    int total_volume() const;
    double adhesion_energy(int i, int new_sigma);
    double volume_energy(int sigma_old, int new_sigma);
    void set_temperature( double temperature );

private:
    std::vector<int> lattice;
    std::vector<int> volume;
    std::vector<int> V0;
    std::vector<double> J;

    double lambdaV;
    double temperature;

    std::mt19937& gen;
    std::uniform_int_distribution<int> dist_site;
    std::uniform_real_distribution<double> dist_prob;
    int idx(int x, int y) const;
    int random_neighbor(int x, int y);
};
