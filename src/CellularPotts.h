#include <vector>
#include <random>
#include <cmath>
#include <iostream>

class CellularPotts {
public:
    CellularPotts(int L, int num_cells, std::mt19937 gen_)
        : L(L), lattice(L*L), 
          volume(num_cells, 0),
          V0(num_cells, 100),
          lambdaV(50.0),
          temperature(10.0),
          gen(gen_),
          dist_site(0, L-1),
          dist_prob(0.0, 1.0)
    {
        // Initialize lattice with random cell IDs
        std::uniform_int_distribution<int> dist_cell(0, num_cells-1);
        for (int i = 0; i < L*L; ++i) {
            lattice[i] = dist_cell(gen);
            volume[lattice[i]]++;
        }

        // Example adhesion matrix
        J.resize(num_cells, std::vector<double>(num_cells, 16.0));
    }

    void monte_carlo_step() {
        for (int k = 0; k < L*L; ++k)
            attempt_copy();
    }

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

    void attempt_copy() ;
};
