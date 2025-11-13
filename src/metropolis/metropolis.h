//
// Created by ozdalkiran-l on 11/13/25.
//

#ifndef LQCD_ECMC_EMBEDDING_METROPOLIS_H
#define LQCD_ECMC_EMBEDDING_METROPOLIS_H

#include "../su3utils/su3utils.h"
#include "../lattice/lattice.h"


void metropolis_sweep(vector<Complex> &links, const Lattice &lat, double beta, mt19937_64 &rng, vector<SU3> &set, size_t &accepted, size_t &proposed, int n_hits);

vector<double> metropolis_samples(vector<Complex> &links, const Lattice &lat, double beta, double epsilon, int n_set, int n_meas, int n_sweeps_meas, int n_hits, int n_burnin, size_t &accepted, size_t &proposed, mt19937_64 &rng);

#endif //LQCD_ECMC_EMBEDDING_METROPOLIS_H