//
// Created by ozdalkiran-l on 11/13/25.
//

#include "observables.h"
#include <algorithm>
#include "../lattice/lattice.h"
#include "../su3utils/su3utils.h"


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

SU3 clover_site(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, int nu) {
    //Computes G_{\mu\nu}(site) clover
    if (mu==nu) cerr << "mu = nu => G(site, mu, nu) = 0" << endl;
    size_t x = site;
    size_t xpmu = lat.neighbor[x][mu][0]; //x+mu
    size_t xpnu = lat.neighbor[x][nu][0]; //x+nu
    size_t xmmu = lat.neighbor[x][mu][1]; //x-mu
    size_t xmnu = lat.neighbor[x][nu][1]; //x-nu
    size_t xmmupnu = lat.neighbor[xmmu][nu][0]; //x-mu+nu
    size_t xmmumnu = lat.neighbor[xmmu][nu][1]; //x-mu-nu
    size_t xpmumnu = lat.neighbor[xpmu][nu][1]; //x+mu-nu
    SU3 clover = SU3::Zero();
    clover += view_link_const(links, x, mu) * view_link_const(links, xpmu, nu) * view_link_const(links, xpnu, mu).adjoint() * view_link_const(links, x, nu).adjoint();
    clover += view_link_const(links,x, nu) * view_link_const(links,xmmupnu, mu).adjoint() * view_link_const(links,xmmu,nu).adjoint() * view_link_const(links,xmmu,mu).adjoint();
    clover += view_link_const(links,xmmu, mu).adjoint() * view_link_const(links,xmmumnu, nu).adjoint() * view_link_const(links,xmmumnu, mu) * view_link_const(links,xmnu, nu);
    clover += view_link_const(links,xmnu, nu).adjoint() * view_link_const(links,xmnu, mu) * view_link_const(links,xpmumnu, nu) * view_link_const(links,x, mu).adjoint();
    clover = 0.25 * clover.imag();
    return clover;
}

double local_topo_charge_clover(const vector<Complex> &links, const Lattice &lat, size_t site) {
    //Computes the local clover topological charge at site
    double q_clover = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho =0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    if (levi_civita(mu,nu,rho,sigma) != 0) {
                        double tr = (clover_site(links, lat, site, mu, nu) * clover_site(links, lat, site, rho, sigma)).trace().real();
                        q_clover += levi_civita(mu, nu, rho, sigma) * tr;
                    }
                }
            }
        }
    }
    q_clover *= 1.0/(32 * M_PI * M_PI);
    return q_clover;
}

double topo_charge_clover(const vector<Complex> &links, const Lattice &lat) {
    double q = 0.0;
    for (size_t site = 0; site < lat.V; site++) {
        q += local_topo_charge_clover(links, lat, site);
    }
    return q;
}

