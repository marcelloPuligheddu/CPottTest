#include <vector>
#include <random>
#include <cmath>
#include <iostream>

CellularPotts :: CellularPotts(int L, int num_cells, std::std::mt19937 gen)
        : L(L), lattice(L*L), 
          volume(num_cells, 0),
          V0(num_cells, 100),
          lambdaV(50.0),
          temperature(10.0),
          gen(gen),
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

    void CellularPotts::monte_carlo_step() {
        for (int k = 0; k < L*L; ++k)
            attempt_copy();
    }

    inline int CellularPotts::idx(int x, int y) const {
        return x + L*y;
    }

    int CellularPotts::neighbor(int x, int y) {
        int dir = gen() % 4;
        if (dir == 0) x = (x+1) % L;
        if (dir == 1) x = (x-1+L) % L;
        if (dir == 2) y = (y+1) % L;
        if (dir == 3) y = (y-1+L) % L;
        return idx(x,y);
    }

    double CellularPotts::adhesion_energy(int i, int new_sigma) {
        int x = i % L;
        int y = i / L;
        double dH = 0.0;

        const int sigma_old = lattice[i];

        // 4-neighborhood
        int neigh[4] = {
            idx((x+1)%L, y),
            idx((x-1+L)%L, y),
            idx(x, (y+1)%L),
            idx(x, (y-1+L)%L)
        };

        for (int n : neigh) {
            int sigma_n = lattice[n];

            dH += J[new_sigma][sigma_n] * (new_sigma != sigma_n);
            dH -= J[sigma_old][sigma_n] * (sigma_old != sigma_n);
        }

        return dH;
    }

    double CellularPotts::volume_energy(int sigma_old, int sigma_new) {
        double dH = 0.0;

        // remove from old
        int V_old = volume[sigma_old];
        dH += lambdaV * ((V_old-1 - V0[sigma_old])*(V_old-1 - V0[sigma_old])
                       - (V_old - V0[sigma_old])*(V_old - V0[sigma_old]));

        // add to new
        int V_new = volume[sigma_new];
        dH += lambdaV * ((V_new+1 - V0[sigma_new])*(V_new+1 - V0[sigma_new])
                       - (V_new - V0[sigma_new])*(V_new - V0[sigma_new]));

        return dH;
    }

    void CellularPotts::attempt_copy() {
        int x = dist_site(gen);
        int y = dist_site(gen);
        int i = idx(x,y);

        int j = neighbor(x,y);

        int sigma_old = lattice[i];
        int sigma_new = lattice[j];

        if (sigma_old == sigma_new) return;

        double dH = 0.0;
        dH += adhesion_energy(i, sigma_new);
        dH += volume_energy(sigma_old, sigma_new);

        if (dH <= 0.0 || dist_prob(gen) < std::exp(-dH/temperature)) {
            lattice[i] = sigma_new;
            volume[sigma_old]--;
            volume[sigma_new]++;
        }
    }
};
