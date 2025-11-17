//
// Created by ozdalkiran-l on 11/13/25.
//
#include "lattice.h"
#include "../su3utils/su3utils.h"

using Complex = std::complex<double>;
using SU3 = Eigen::Matrix3cd;

using namespace std;

Lattice::Lattice(int Nx_, int Ny_, int Nz_, int Nt_) : Nx(Nx_), Ny(Ny_), Nz(Nz_), Nt(Nt_) {
    V = static_cast<size_t>(Nx) * Ny * Nz * Nt;
    neighbor.resize(V);
    //Calcul des voisins de chaque site avec conditions aux bords périodiques
    for (int t = 0; t < Nt; t++) {
        for (int z = 0; z < Nz; z++) {
            for (int y = 0; y < Ny; y++) {
                for (int x = 0; x < Nx; x++) {
                    size_t i = index(x, y, z, t);
                    neighbor[i][0][0] = index((x + 1) % Nx, y, z, t);
                    neighbor[i][0][1] = index((x - 1 + Nx) % Nx, y, z, t);
                    neighbor[i][1][0] = index(x, (y + 1) % Ny, z, t);
                    neighbor[i][1][1] = index(x, (y - 1 + Ny) % Ny, z, t);
                    neighbor[i][2][0] = index(x, y, (z + 1) % Nz, t);
                    neighbor[i][2][1] = index(x, y, (z - 1 + Nz) % Nz, t);
                    neighbor[i][3][0] = index(x, y, z, (t + 1) % Nt);
                    neighbor[i][3][1] = index(x, y, z, (t - 1 + Nt) % Nt);
                }
            }
        }
    }
    staples.resize(V);
    //Calcul des liens de chaque staple autour de chaque lien
    for (int t = 0; t < Nt; t++) {
        for (int z = 0; z < Nz; z++) {
            for (int y = 0; y < Ny; y++) {
                for (int x = 0; x < Nx; x++) {
                    size_t site = index(x, y, z, t); //x
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;

                            size_t xmu = neighbor[site][mu][0]; //x+mu
                            size_t xnu = neighbor[site][nu][0]; //x+nu
                            size_t xmunu = neighbor[xmu][nu][1]; //x+mu-nu
                            size_t xmnu = neighbor[site][nu][1]; //x-nu

                            staples[site][mu][j][0] = {xmu, nu};
                            staples[site][mu][j][1] = {xnu, mu};
                            staples[site][mu][j][2] = {site, nu};

                            staples[site][mu][j + 1][0] = {xmunu, nu};
                            staples[site][mu][j + 1][1] = {xmnu, mu};
                            staples[site][mu][j + 1][2] = {xmnu, nu};

                            j += 2;
                        }
                    }
                }
            }
        }
    }
}

void hot_start(vector<Complex> &links, const Lattice &lat, mt19937_64 &rng) {
    //Effectue un hot start sur une configuration de jauge
    //std::random_device rd;
    //std::mt19937_64 rng(rd());
    for (size_t site = 0; site < lat.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            view_link(links, site, mu) = random_su3(rng);
        }
    }
}

void cold_start(vector<Complex> &links, const Lattice &lat) {
    //Effectue un cold start sur une configuration de jauge
    for (size_t site = 0; site < lat.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            view_link(links, site, mu) = Eigen::Matrix3cd::Identity();
        }
    }
}

void compute_staple(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, SU3 &staple) {
    //Calcule la somme des staples sur un lien de jauge
    staple.setZero();
    for (int nu = 0; nu < 4; nu++) {
        if (nu == mu) {
            continue;
        }
        size_t x = site; //x
        size_t xmu = lat.neighbor[x][mu][0]; //x+mu
        size_t xnu = lat.neighbor[x][nu][0]; //x+nu
        size_t xmunu = lat.neighbor[xmu][nu][1]; //x+mu-nu
        size_t xmnu = lat.neighbor[x][nu][1]; //x-nu
        auto U0 = view_link_const(links, xmu, nu);
        auto U1 = view_link_const(links, xnu, mu);
        auto U2 = view_link_const(links, x, nu);
        staple += U0 * U1.adjoint() * U2.adjoint();
        auto V0 = view_link_const(links, xmunu, nu);
        auto V1 = view_link_const(links, xmnu, mu);
        auto V2 = view_link_const(links, xmnu, nu);
        staple += V0.adjoint() * V1.adjoint() * V2;
    }
}