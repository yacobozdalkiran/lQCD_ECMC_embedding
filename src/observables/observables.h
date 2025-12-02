//
// Created by ozdalkiran-l on 11/13/25.
//

#ifndef LQCD_ECMC_EMBEDDING_OBSERVABLES_H
#define LQCD_ECMC_EMBEDDING_OBSERVABLES_H

#include "../lattice/lattice.h"

double mean_plaquette(const vector<Complex> &links, const Lattice &lat);

struct PlaquetteStats {
    double mean;
    double stddev;
};

PlaquetteStats plaquette_stats(const vector<Complex> &links, const Lattice &lat);

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

//Fonctions pour check amélioration refresh R
//qtr = quasi trace (i.e. p00.imag() +p11.imag())
//Remplacement : qtr = (p11.imag() - p00.imag()) / (p00.imag() + p11.imag())

inline double mean(const std::array<double,4>& a) {
    double sum = std::accumulate(a.begin(), a.end(), 0.0);
    return sum / a.size();   // ici size() = 4
}

inline double variance(const std::array<double,4>& a) {
    double m = mean(a);

    double sum_sq = 0.0;
    for (double x : a)
        sum_sq += (x - m) * (x - m);

    return sum_sq / a.size();   // variance population
}

void measure_qtr_one_link_one_j(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, int j, vector<SU3> &P,array<double, 4> &qtr_link_j,array<double, 4> &stat_qtr_link_j);

void measure_qtr(const vector<Complex> &links, const Lattice &lat, vector<SU3> &P, vector<double> &mins, vector<double> &maxs, vector<double> &means, vector<double> &vars);

#endif //LQCD_ECMC_EMBEDDING_OBSERVABLES_H