//
// Created by ozdalkiran-l on 11/13/25.
//

#ifndef LQCD_ECMC_EMBEDDING_OBSERVABLES_H
#define LQCD_ECMC_EMBEDDING_OBSERVABLES_H

#include "../lattice/lattice.h"

//Plaquette moyenne
double mean_plaquette(const vector<Complex> &links, const Lattice &lat);
struct PlaquetteStats {
    double mean;
    double stddev;
};
PlaquetteStats plaquette_stats(const vector<Complex> &links, const Lattice &lat);

//Action Wilson
double wilson_action(const vector<Complex> &links, const Lattice &lat, double beta);

//Charge topologique
SU3 clover_site(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, int nu);
inline int levi_civita(int mu, int nu, int rho, int sigma) {
    if (mu == nu || mu == rho || mu == sigma || nu == rho || nu == sigma || rho == sigma) return 0;
    int inv = 0;
    inv += (mu > nu);
    inv += (mu > rho);
    inv += (mu > sigma);
    inv += (nu > rho);
    inv += (nu > sigma);
    inv += (rho > sigma);

    return (inv & 1) ? -1 : +1; //Parity test on inv
}
double local_topo_charge_clover(const vector<Complex> &links, const Lattice &lat, size_t site);
double topo_charge_clover(const vector<Complex> &links, const Lattice &lat);

#endif //LQCD_ECMC_EMBEDDING_OBSERVABLES_H