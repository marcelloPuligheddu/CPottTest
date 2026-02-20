#include <vector>
#include <random>
#include <cmath>
#include <iostream>


#include"CellularPotts.h"

template<int L, int K >
double CellularPotts<L,K>::CellularPotts( std::mt19937& gen)
        : lattice(L*L), 
          volume(K, 0),
          V0(K, 100),
          J(K*K, 100),
          lambdaV(50.0),
          temperature(10.0),
          gen(gen),
          dist_site(0, L-1),
          dist_prob(0.0, 1.0) {
        // Initialize lattice with random cell IDs
        std::uniform_int_distribution<int> dist_cell(0, K-1);
        for (int i = 0; i < L*L; ++i) {
            lattice[i] = dist_cell(gen);
            volume[lattice[i]]++;
        }

}

template<int L, int K >
inline int CellularPotts<L,K>::idx(int x, int y) const { return x + L*y; }

template<int L, int K >
int CellularPotts<L,K>::random_neighbor(int x, int y) {
        int dir = gen() % 4;
        if (dir == 0) x = (x+1) % L;
        if (dir == 1) x = (x-1+L) % L;
        if (dir == 2) y = (y+1) % L;
        if (dir == 3) y = (y-1+L) % L;
        return idx(x,y);
}


template<int L, int K >
void CellularPotts<L,K>::monte_carlo_step() { for (int k = 0; k < L*L; ++k) attempt_copy(); }

template<int L, int K >
double CellularPotts<L,K>::compute_delta(int i, int j) const {
    int sigma_old = lattice[i];
    int sigma_new = lattice[j];

    if (sigma_old == sigma_new) return 0.0;

    double dH = 0.0;

    // adhesion
    int x = i % L;
    int y = i / L;

    int neigh[4] = {
        idx((x+1)%L, y),
        idx((x-1+L)%L, y),
        idx(x, (y+1)%L),
        idx(x, (y-1+L)%L)
    };

    for (int n : neigh) {
        int sigma_n = lattice[n];
        dH += J[sigma_new*K+sigma_n] * (sigma_new != sigma_n);
        dH -= J[sigma_old*K+sigma_n] * (sigma_old != sigma_n);
    }

    // volume
    int V_old = volume[sigma_old];
    int V_new = volume[sigma_new];

    dH += lambdaV * (
        (V_old-1 - V0[sigma_old])*(V_old-1 - V0[sigma_old])
      - (V_old   - V0[sigma_old])*(V_old   - V0[sigma_old])
    );

    dH += lambdaV * (
        (V_new+1 - V0[sigma_new])*(V_new+1 - V0[sigma_new])
      - (V_new   - V0[sigma_new])*(V_new   - V0[sigma_new])
    );

    return dH;
}

template<int L, int K >
void CellularPotts<L,K>::apply_move(int i, int j)
{
    int sigma_old = lattice[i];
    int sigma_new = lattice[j];

    lattice[i] = sigma_new;
    volume[sigma_old]--;
    volume[sigma_new]++;
}

template<int L, int K >
double CellularPotts<L,K>::adhesion_energy(int i, int new_sigma) {
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

            dH += J[new_sigma*K+sigma_n] * (new_sigma != sigma_n);
            dH -= J[sigma_old*K+sigma_n] * (sigma_old != sigma_n);
        }

        return dH;
}

template<int L, int K >
double CellularPotts<L,K>::volume_energy(int sigma_old, int sigma_new) {
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

template<int L, int K >
bool CellularPotts<L,K>::accept(double dH)
{
    if (dH <= 0.0)
        return true;

    return dist_prob(gen) < std::exp(-dH / temperature);
}

template<int L, int K >
bool CellularPotts<L,K>::attempt_copy()
{
    int x = dist_site(gen);
    int y = dist_site(gen);

    int i = idx(x,y);
    int j = random_neighbor(x,y);

    double dH = compute_delta(i, j);

    if (accept(dH)) {
        apply_move(i, j);
        return true;
    }

    return false;
}

template<int L, int K >
double CellularPotts<L,K>::total_energy() const {
        double H = 0.0;
    
        // --- Adhesion term ---
        for (int y = 0; y < L; ++y)
        {
            for (int x = 0; x < L; ++x)
            {
                int i = idx(x, y);
                int sigma_i = lattice[i];
    
                // Right neighbor
                int xr = (x + 1) % L;
                int ir = idx(xr, y);
                int sigma_r = lattice[ir];
    
                if (sigma_i != sigma_r)
                    H += J[sigma_i][sigma_r];
    
                // Up neighbor
                int yu = (y + 1) % L;
                int iu = idx(x, yu);
                int sigma_u = lattice[iu];
    
                if (sigma_i != sigma_u)
                    H += J[sigma_i][sigma_u];
            }
        }
    
        // --- Volume constraint ---
        for (size_t s = 0; s < volume.size(); ++s)
        {
            double dv = volume[s] - V0[s];
            H += lambdaV * dv * dv;
        }
    
        return H;
}

template<int L, int K >
void CellularPotts<L,K>::set_temperature( double new_temperature ){ temperature = new_temperature ;  }

template<int L, int K >
int CellularPotts<L,K>::total_volume() const {
        int sum = 0;
        for (int v : volume)
            sum += v;
        return sum;
}

