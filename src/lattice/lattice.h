//
// Created by ozdalkiran-l on 11/13/25.
//

#ifndef LQCD_ECMC_EMBEDDING_LATTICE_H
#define LQCD_ECMC_EMBEDDING_LATTICE_H

#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <complex>
#include <random>
#include <string>

using Complex = std::complex<double>;
using SU3 = Eigen::Matrix3cd;

using namespace std;

struct Lattice {
    int Nx, Ny, Nz, Nt;
    size_t V;
    vector<array<array<size_t,2>,4>> neighbor;
    vector<array<array<array<pair<size_t,int>,3>,6>,4>> staples;

    Lattice(int Nx_, int Ny_, int Nz_, int Nt_);
    inline size_t index(int x, int y, int z, int t) const {
        return ((size_t(t)*Nz + z) * Ny + y)*Nx + x;
    }
};

inline Eigen::Map<SU3> view_link(vector<Complex> &links, size_t site, int mu) {
    //Permet de mapper les éléments d'un vecteur de complex en matrices SU3 Eigen (permet d'utiliser Eigen avec contiguité mémoire)
    return Eigen::Map<SU3>(&links[(site * 4 + mu) * 9]);
}

inline Eigen::Map<const SU3> view_link_const(const vector<Complex> &links, size_t site, int mu) {
    //Pareil que view_link mais en const pour safety
    return Eigen::Map<const SU3>(&links[(site * 4 + mu) * 9]);
}

void hot_start(vector<Complex> &links, const Lattice &lat, mt19937_64 &rng);

void cold_start(vector<Complex> &links, const Lattice &lat);

void compute_staple(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, SU3 &staple);

#endif //LQCD_ECMC_EMBEDDING_LATTICE_H