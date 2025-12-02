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
    double qclover = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho =0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    if (levi_civita(mu,nu,rho,sigma) != 0) {
                        double tr = (clover_site(links, lat, site, mu, nu) * clover_site(links, lat, site, rho, sigma)).trace().real();
                        qclover += levi_civita(mu, nu, rho, sigma) * tr;
                    }
                }
            }
        }
    }
    qclover *= 1.0/(32 * M_PI * M_PI);
    return qclover;
}

double topo_charge_clover(const vector<Complex> &links, const Lattice &lat) {
    double q = 0.0;
    for (size_t site = 0; site < lat.V; site++) {
        q += local_topo_charge_clover(links, lat, site);
    }
    return q;
}

void measure_qtr_one_link_one_j(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, int j, vector<SU3> &P,array<double, 4> &qtr_link_j,array<double, 4> &stat_qtr_link_j) {
    SU3 U0 = view_link_const(links, site, mu);
    auto links_staple_j = lat.staples[site][mu][j]; //Liste des 3 liens de la staple j
    SU3 U1 = view_link_const(links, links_staple_j[0].first, links_staple_j[0].second); //Les 3 matrices SU3 associées
    SU3 U2 = view_link_const(links, links_staple_j[1].first, links_staple_j[1].second);
    SU3 U3 = view_link_const(links, links_staple_j[2].first, links_staple_j[2].second);

    array<pair<size_t, int>,4> links_plaquette_j; //On rajoute le lien actuel
    links_plaquette_j[0] = make_pair(site, mu);
    links_plaquette_j[1] = links_staple_j[0];
    links_plaquette_j[2] = links_staple_j[1];
    links_plaquette_j[3] = links_staple_j[2];

    if (j%2 == 0) { //Forward plaquette
        P[0] = U0 * U1 * U2.adjoint() * U3.adjoint();
        P[1] = U1 * U2.adjoint() * U3.adjoint() * U0;
        P[2] = U2 * U1.adjoint() * U0.adjoint()* U3;
        P[3] = U3 * U2 * U1.adjoint() * U0.adjoint();
    }
    else { //Backward plaquette
        P[0] = U0 * U1.adjoint() * U2.adjoint() * U3;
        P[1] = U1 * U0.adjoint() * U3.adjoint() * U2;
        P[2] = U2 * U1 * U0.adjoint()* U3.adjoint();
        P[3] = U3 * U0 * U1.adjoint() * U2.adjoint();
    }

    qtr_link_j[0] = (-P[0](0,0).imag() + P[0](1,1).imag())/(P[0](0,0).imag() + P[0](1,1).imag());
    qtr_link_j[1] = (-P[1](0,0).imag() + P[1](1,1).imag())/(P[1](0,0).imag() + P[1](1,1).imag());
    qtr_link_j[2] = (-P[2](0,0).imag() + P[2](1,1).imag())/(P[2](0,0).imag() + P[2](1,1).imag());
    qtr_link_j[3] = (-P[3](0,0).imag() + P[3](1,1).imag())/(P[3](0,0).imag() + P[3](1,1).imag());
    auto [mn, mx] = std::minmax_element(qtr_link_j.begin(), qtr_link_j.end());
    stat_qtr_link_j[0] = *mn;
    stat_qtr_link_j[1] = *mx;
    stat_qtr_link_j[2] = mean(qtr_link_j);
    stat_qtr_link_j[3] = variance(qtr_link_j);
}

void measure_qtr(const vector<Complex> &links, const Lattice &lat, vector<SU3> &P, vector<double> &mins, vector<double> &maxs, vector<double> &means, vector<double> &vars){
    size_t N_stats = lat.V * 4 * 6;
    mins.resize(N_stats);
    maxs.resize(N_stats);
    means.resize(N_stats);
    vars.resize(N_stats);
    array<double, 4> qtr_link_j;
    array<double, 4> stat_qtr_link_j;
    for (size_t site = 0; site < lat.V; site++) {
        for (int mu = 0; mu<4; mu++) {
            for (int j = 0; j < 6; j++) {
                measure_qtr_one_link_one_j(links, lat, site, mu, j, P, qtr_link_j,stat_qtr_link_j);
                size_t index = (site * 4 * 6) + mu * 6 + j;
                mins[index] = stat_qtr_link_j[0];
                maxs[index] = stat_qtr_link_j[1];
                means[index] = stat_qtr_link_j[2];
                vars[index] = stat_qtr_link_j[3];
            }
        }
    }
}

