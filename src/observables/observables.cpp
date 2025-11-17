//
// Created by ozdalkiran-l on 11/13/25.
//

#include "observables.h"

double mean_plaquette(const vector<Complex> &links, const Lattice &lat) {
    //Calcule la plaquette moyenne d'une configuration
    double result = 0;
    size_t counter = 0;
    SU3 staple;
    for (size_t site = 0; site < lat.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            compute_staple(links, lat, site, mu, staple);
            result += (view_link_const(links, site, mu)*staple).trace().real()/3.0;
            counter += 6; //6 plaquettes par site
        }
    }
    return result / static_cast<double>(counter);
}

PlaquetteStats plaquette_stats(const vector<Complex> &links, const Lattice &lat) {
    //Calcule la plaquette moyenne et la variance sur la plaquette moyenne au sein d'une configuration
    double sum = 0.0;
    double sum2 = 0.0;
    size_t counter = 0;
    SU3 staple;

    for (size_t site = 0; site < lat.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            compute_staple(links, lat, site, mu, staple);
            double p = (view_link_const(links, site, mu) * staple).trace().real() / 18.0; //6 plaquettes par site, facteur 3 car SU(3)
            sum += p;
            sum2 += p * p;
            counter += 1;
        }
    }

    double mean = sum / static_cast<double>(counter);
    double var = (sum2 / static_cast<double>(counter)) - mean * mean;
    double stddev = std::sqrt(var / static_cast<double>(counter)); // erreur sur la moyenne (σ/√N)

    return {mean, stddev};
}
