//
// Created by ozdalkiran-l on 11/13/25.
//

#ifndef LQCD_ECMC_EMBEDDING_ECMC_H
#define LQCD_ECMC_EMBEDDING_ECMC_H

#include "../lattice/lattice.h"

void compute_list_staples(const vector<Complex> &links, const Lattice &lat, size_t site, int mu,
    array<SU3,6> &list_staple);

void compute_reject(double A, double B, double &gamma, double &reject, int epsilon);

void compute_reject_angles(const vector<Complex> &links, size_t site, int mu, const array<SU3,6> &list_staple,
    const SU3 &R, int epsilon, const double &beta, array<double,6> &reject_angles, mt19937_64 &rng);

int selectVariable(const vector<double> &probas, mt19937_64 &rng);

pair<pair<size_t, int>,int> lift(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, int j,
    const SU3 &R, const SU3 &lambda_3, mt19937_64 &rng);

pair<pair<size_t, int>,int> lift_improved(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, int j,
    SU3 &R, const SU3 &lambda_3, mt19937_64 &rng, const vector<SU3> &set);

void ecmc_update(vector<Complex> &links, size_t site, int mu, double theta, int epsilon, const SU3 &R);

vector<double> ecmc_samples(vector<Complex> &links, const Lattice &lat, double beta, int N_samples,
    double param_theta_sample, double param_theta_refresh, mt19937_64 &rng, bool poisson);

vector<double> ecmc_samples_improved(vector<Complex> &links, const Lattice &lat, double beta, int N_samples,
    double param_theta_sample, double param_theta_refresh, mt19937_64 &rng, bool poisson, double epsilon_set);

//ECMC temps discret
int update_until_reject_d(vector<Complex> &links, size_t site, int mu, const array<SU3, 6> &list_staple,
    const SU3 &R, int epsilon, const double &beta, const double &eta, mt19937_64 &rng);
#endif //LQCD_ECMC_EMBEDDING_ECMC_H