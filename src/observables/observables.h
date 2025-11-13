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


#endif //LQCD_ECMC_EMBEDDING_OBSERVABLES_H